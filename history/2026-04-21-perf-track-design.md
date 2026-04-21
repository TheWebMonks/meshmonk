# Design: Performance Track — Benchmark Harness + Core Parallelization

**Date:** 2026-04-21
**Status:** Draft
**Decisions:** [ADR-002](../docs/decisions/ADR-002-perf-track.md)
**Origin:** Post-e2e-test triage scan (2026-04-21). Four related beads surfaced — `meshmonk-modernization-44m` (benchmark harness), `-9f5` (KD-tree parallelization), `-6a5` (cache temporaries), `-bdt` (smoothing parallelization) — that should move as a single sequenced epic rather than four independent changes.
**Supersedes:** —
**Parent:** [ADR-001 modernization](../docs/decisions/ADR-001-meshmonk-modernization.md)

---

## Problem

MeshMonk's Python layer is stable and green on CI, but the registration hot paths have never been measured or parallelized. Profiling by inspection identifies three candidate bottlenecks:

1. **KD-tree neighbour queries** (`library/src/NeighbourFinder.hpp:157-181`) — sequential loop over ~50K query points, called every correspondence pass. This is the dominant cost of rigid + nonrigid iteration.
2. **Per-iteration temporary allocations** in nonrigid/viscoelastic updates (`NonrigidRegistration.cpp`, `ViscoElasticTransformer.cpp:141-248`, `InlierDetector.cpp:84-127`) — buffers like `regularizedForceField`, `temporaryDisplacementField`, `tempInlierWeights` are allocated fresh every ICP iteration.
3. **Smoothing loops** (`ViscoElasticTransformer.cpp:166-239`, `InlierDetector.cpp:90-121`) — per-vertex independent weighted averages that trivially parallelize.

The deeper issue is that we have **no benchmark harness**, so any perf claim today would be guessed. Optimizing without measurement is how the codebase accumulates untested "speedups" that drift or regress. The harness is a strict prerequisite.

OpenMP is not currently linked (`CMakeLists.txt:44-45` has a comment confirming `vci_openmp()` is never invoked). Adding OpenMP is a single `find_package(OpenMP)` + `target_link_libraries(... OpenMP::OpenMP_CXX)` — the toolchain shift is minor but it IS a toolchain shift and belongs in the ADR.

## Goal

Deliver measured, regression-guarded performance improvements to the three hot paths above, in a sequence that keeps every merge safe:

- Benchmark harness lands first and captures baseline numbers.
- Each optimization merges only after demonstrating a measured wall-clock win on the harness AND preserving the existing e2e test numerics (within tolerance — see D3 on parallel reduction non-determinism).
- No change to the public API, no algorithmic changes, no new dependencies except OpenMP.

## Non-Goals

- C++ unit-test harness (tracked separately in `meshmonk-modernization-0bs`).
- API/docstring/type-hint polish (tracked separately in `-hnf`, `-xvu`, `-h31`).
- Algorithmic changes (e.g. switching KD-tree backend, different smoothing scheme).
- GPU/SIMD/AVX-512 — out of scope; OpenMP CPU threading only.
- Thread-pool tuning knobs on the Python surface. `OMP_NUM_THREADS` env var is sufficient for v0.x.
- Multi-process parallelism (`multiprocessing`, joblib) — we want intra-call parallelism, not batch.

---

## Proposed Architecture

The track is one epic with four sequenced children:

```
perf-track (epic)
├── 44m  Benchmark harness (pytest-benchmark)       [P1, blocks all others]
├── 9f5  Parallelize KD-tree neighbour queries       [P1, depends on 44m]
├── 6a5  Cache per-iteration temporaries              [P2, depends on 44m]
└── bdt  Parallelize smoothing loops                 [P2, depends on 44m, 9f5]
```

`bdt` sequences after `9f5` so we only add OpenMP CMake wiring once. `6a5` is independent of the other two (no OpenMP, pure allocation hoisting) and can interleave.

### Component: Benchmark harness (`44m`)

`benchmarks/` top-level directory using `pytest-benchmark`. One file per registration type (`bench_rigid.py`, `bench_nonrigid.py`, `bench_pyramid.py`) with a parametrized mesh-size fixture covering (1K, 5K, 50K) vertices. Inputs synthesized from `data/faces/Template.obj` via downsampling (deterministic) + a known rigid perturbation.

`pytest-benchmark` is **not** currently in the dev dependency set — `44m` adds it to `[project.optional-dependencies].dev` in `pyproject.toml`. It is NOT added to cibuildwheel `test-requires`; wheel CI does not run benchmarks. The regression gate is a separate nightly workflow.

Harness runs two modes:
- **Local**: full `pytest-benchmark compare` workflow; writes `.benchmarks/` results to user's machine (gitignored).
- **CI regression gate**: a single nightly workflow runs a curated subset (rigid @ 5K, nonrigid @ 5K, pyramid @ 5K) and compares against a committed `benchmarks/baseline.json`. Fails if rigid/nonrigid scenarios regress by >20% wall-clock; pyramid uses a 30% tolerance because its iteration count varies with annealing convergence and is noisier on CI runners. `pytest-benchmark` is configured with `rounds >= 5` and median-of-rounds to dampen outliers.

**Baseline isolation:** baseline runs pin `OMP_NUM_THREADS=1` and `EIGEN_DONT_PARALLELIZE=1` so the baseline captures strictly single-threaded MeshMonk behavior. Post-change benchmarks run with whatever is native, so the reported speedup is the *incremental* effect of MeshMonk-level OpenMP pragmas — not contaminated by Eigen's own opportunistic threading.

Baseline JSON is refreshed explicitly (PR-gated), not auto-updated on green builds. This keeps the baseline meaningful as a regression anchor rather than drifting with noise.

### Component: KD-tree parallelization (`9f5`)

`NeighbourFinder<VecMatType>::update()` becomes `#pragma omp parallel for` over the query loop at `NeighbourFinder.hpp:157`. Requires:

- Making `queriedFeature`, `neighbourIndices`, `neighbourSquaredDistances`, and the `knnResultSet` per-thread-local (declare inside the loop body, not outside).
- Verifying `nanoflann`'s `findNeighbors` is thread-safe on a shared `_kdTree`. nanoflann's `KDTreeSingleIndexAdaptor` is read-only-safe **iff the index is not being rebuilt**. The kd-tree build in `set_source_points` (invoked by `ViscoElasticTransformer::_update_neighbours` on each outer ICP iteration) stays sequential and happens OUTSIDE the parallel region — only the per-query loop is parallelized. `SearchParams(32, 0.0001, true)` with `sorted=true` is safe under parallel queries (sorting is per-query result set, not global).
- Output arrays `_outNeighbourIndices` and `_outNeighbourSquaredDistances` are written at distinct row indices `i` per thread — no race.

**CMake detail:** MSVC's default `/openmp` flag is OpenMP 2.0 and rejects `size_t` loop counters. `library/CMakeLists.txt` must force `/openmp:llvm` on MSVC: `if(MSVC) target_compile_options(meshmonk_lib PRIVATE /openmp:llvm) endif()`. This is a required default, not a fallback.

Expected: 4–8× speedup on correspondence phase on 8-core machines. Verify against baseline before declaring done.

### Component: Cache temporaries (`6a5`)

Hoist repeated allocations out of ICP iteration loops and reuse via `setZero()` / `conservativeResize()`:

- `NonrigidRegistration.cpp:59-61` — `correspondingFeatures`, etc. currently allocated per iteration.
- `ViscoElasticTransformer.cpp` — `regularizedForceField`, `unregulatedDisplacementField`, `temporaryDisplacementField` (lines 162, 205, 248).
- `InlierDetector.cpp:86` — `tempInlierWeights`.

Move these to member fields on the respective classes, sized in `set_input`/`set_output` (matching the existing pattern in `ViscoElasticTransformer::set_output:25-27`, which already resizes `_displacementField`/`_oldDisplacementField` on each set), and zeroed at top of each iteration. **Do NOT** size-on-first-call: that pattern breaks if a registration instance is reused with a different-sized floating mesh (e.g. pyramid layers with varying vertex counts). Today `pyramid_register` instantiates a fresh registration per layer, but the buffer sizing must not silently depend on that — size-in-setter is the robust pattern. A small unit test exercising back-to-back `update()` calls with different input sizes catches regressions.

Expected: 5–15% wall-clock on pyramid (which does the most iterations).

No API change. Numeric output should be bitwise identical (same operations, different memory).

### Component: Smoothing parallelization (`bdt`)

Three loops, all with independent per-vertex work (`#pragma omp parallel for schedule(static)` with chunk size ≥ 16 rows to mitigate false-sharing on contiguous `Vec3Mat` member buffers — each row is 12 bytes, 16 rows ≈ 192 bytes clears a 64-byte cache line comfortably):
- `ViscoElasticTransformer::_update_viscously` smoothing loop (`:167`)
- `ViscoElasticTransformer::_update_elastically` inner loop (`:214`)
- `InlierDetector::_smooth_inlier_weights` inner loop (`:95`)

Each loop reads from one buffer and writes to a different buffer at row `i` — no aliasing. **Invariant (code-review-enforced):** the read-source buffer and the write-destination buffer under 6a5+bdt remain distinct matrices. The existing code is already careful about this (`forceField` → `regularizedForceField`, `_displacementField` → `unregulatedDisplacementField`); after hoisting to member fields in 6a5, both endpoints remain distinct — maintain that. Thread-local scalars (`vectorAverage`, `sumWeights`) stay loop-local by declaration inside the body.

The outlier-diffusion loop (`ViscoElasticTransformer::_update_outlier_transformation`, `:258`) is structurally parallelizable but is deliberately **excluded from `bdt`** because it runs fewer iterations (`_outlierDiffusionIterations`) and has an inlier-skip branch that reduces per-row work — likely below the threading-overhead threshold. Reconsider if profiling shows it hot.

**Parallel-vs-sequential parity test:** `9f5` and `bdt` each land with a test that runs the affected registration path under `OMP_NUM_THREADS=1` vs default threads and asserts output equivalence within a tight ε (e.g. 1e-5 on aligned vertex coordinates). The existing e2e tolerances are too coarse to detect subtle loop-variable races — this test fills the gap.

Expected: 2–4× on smoothing phases.

---

## Key Design Decisions

See ADR-002 for full rationale.

- **D1**: Harness-first sequencing — no optimization merges without measured improvement on the harness.
- **D2**: OpenMP over `std::execution` / TBB / thread-pool library.
- **D3**: Tolerance-based numeric invariance check, not bitwise — parallel float reductions are non-associative.
- **D4**: Per-optimization merge gate — measured speedup threshold + e2e tests green.
- **D5**: `OMP_NUM_THREADS` is the only public thread knob for v0.x.

---

## Consequences

**Easier:**
- Future perf work has a harness to measure against.
- Regression gate catches accidental slowdowns in unrelated changes.
- KD-tree and smoothing phases scale with core count.

**Harder:**
- CMake build now depends on an OpenMP-capable compiler. GCC/Clang ship it; Apple Clang on macOS needs `libomp` via Homebrew. cibuildwheel already supports this (libomp is brought in via brew in the macOS image); verify before merging `9f5`.
- Numeric outputs on the e2e tests may shift by floating-point epsilon on parallel builds. Existing e2e tolerances should already accommodate this, but we check per-bead.

**Tradeoffs:**
- We accept non-deterministic parallel-reduction output (bitwise) in exchange for speedup. Determinism-sensitive users can set `OMP_NUM_THREADS=1`.
- Baseline JSON as a regression gate ties CI to a specific runner profile. Documented in `benchmarks/README.md`.

## Known Limitations

- The harness measures wall-clock, not throughput or per-stage breakdown. If a future optimization needs per-stage profiling, that's a separate design.
- `50K`-vertex scenario uses the full `Template.obj`. We don't have a >100K-vertex benchmark yet; revisit if large-mesh perf becomes a user complaint.
- OpenMP on Apple Silicon macOS with Accelerate-linked Eigen can interact oddly with Grand Central Dispatch. Benchmark on macOS before declaring `9f5` done.
- Calling meshmonk register functions from multiple Python threads concurrently (e.g. a `ThreadPoolExecutor` batching registrations) oversubscribes OpenMP — each Python worker spawns `N_CPU` OMP threads. Not a correctness bug, but a perf footgun. Users running batch workloads should set `OMP_NUM_THREADS=1` and parallelize at the Python level instead. Documented in `benchmarks/README.md` and `docs/gotchas.md` (once `gw1` lands).
- `benchmarks/baseline.json` must be refreshed on GitHub Actions runner-image upgrades (`ubuntu-latest` rolls roughly yearly) OR every minor release, whichever comes first. Stale baselines turn the regression gate into noise. Track via a scheduled reminder bead.

---

## Open Questions

1. **Baseline refresh cadence** — per minor release, per runner-image upgrade, or both? Known Limitations commits to "whichever comes first"; revisit after a week of CI noise data from 44m.
2. **Should `6a5` land before or after `9f5`?** Both depend only on `44m`. `6a5` is riskier in terms of subtle state bugs (hoisted buffers might leak state across calls on re-entry) but lower in raw speedup. Recommend `9f5` first — bigger win, simpler to review.

