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

**Re-triage after `44m`:** once the baseline lands, review the per-phase breakdown BEFORE starting `9f5`/`6a5`/`bdt`. If KD-tree correspondence isn't the dominant cost we assumed, `9f5`'s 2× threshold may be unachievable even with a correct implementation, and priorities need re-examination. This is a 15-minute gate, not a redesign — just a sanity check that the optimization ordering still makes sense given real numbers.

### Component: Benchmark harness (`44m`)

`benchmarks/` top-level directory using `pytest-benchmark`. One file per registration type (`bench_rigid.py`, `bench_nonrigid.py`, `bench_pyramid.py`) with a parametrized mesh-size fixture covering (1K, 5K, 50K) vertices. Inputs synthesized from `data/Template.obj` via downsampling (deterministic) + a known rigid perturbation.

**Mesh-size availability check:** `44m` must verify `data/Template.obj` has ≥50K vertices before committing to the 50K parameterization. If it falls short, either (a) reduce top end to match (e.g. 1K, 5K, 30K), or (b) introduce a deterministic upsampling step (edge-split subdivision) — NOT a random duplication, since benchmark inputs must be reproducible. The choice is implementation-time, not design-time.

`pytest-benchmark` is **not** currently in the dev dependency set — `44m` adds it to `[project.optional-dependencies].dev` in `pyproject.toml`. It is NOT added to cibuildwheel `test-requires`; wheel CI does not run benchmarks. The regression gate is a separate nightly workflow.

Harness runs two modes:
- **Local**: full `pytest-benchmark compare` workflow; writes `.benchmarks/` results to user's machine (gitignored).
- **CI regression gate**: a single nightly workflow (`.github/workflows/benchmarks-nightly.yml`) runs a curated subset (rigid @ 5K, nonrigid @ 5K, pyramid @ 5K) and compares against a committed `benchmarks/baseline.json`. Fails if rigid/nonrigid scenarios regress by >20% wall-clock; pyramid uses a 30% tolerance because its iteration count varies with annealing convergence and is noisier on CI runners. `pytest-benchmark` is configured with `rounds >= 5` and median-of-rounds to dampen outliers.

**CI workflow skeleton** (implementation reference, not frozen spec):

```yaml
# .github/workflows/benchmarks-nightly.yml
on:
  schedule: [{cron: "0 6 * * *"}]      # nightly regression gate
  workflow_dispatch:                     # manual refresh / ad-hoc
    inputs:
      refresh_baseline: {type: boolean, default: false}
jobs:
  bench:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - run: pip install -e ".[dev]"
      - name: Run benchmarks
        run: |
          pytest benchmarks/ \
            --benchmark-only \
            --benchmark-json=new.json \
            --benchmark-compare=benchmarks/baseline.json \
            --benchmark-compare-fail=median:20%
      - if: inputs.refresh_baseline
        run: |
          cp new.json benchmarks/baseline.json
          # open PR with the updated baseline (gh pr create)
```

Regression gate and baseline refresh share one workflow — `workflow_dispatch` with `refresh_baseline=true` runs the same steps and opens a PR updating `baseline.json`. This keeps the refresh mechanism defined (not "whoever remembers to paste JSON") while preserving PR-gated review per ADR D6.

**Baseline isolation:** baseline runs pin `OMP_NUM_THREADS=1` and call `Eigen::setNbThreads(1)` via a `benchmarks/conftest.py` fixture, so the baseline captures strictly single-threaded MeshMonk behavior. (`EIGEN_DONT_PARALLELIZE` is a compile-time macro and cannot be toggled at runtime — the env-var form is a no-op.) Post-change benchmarks run with `Eigen::setNbThreads(0)` (default) and whatever `OMP_NUM_THREADS` is native, so the reported speedup is the *incremental* effect of MeshMonk-level OpenMP pragmas — not contaminated by Eigen's own opportunistic threading. In practice, once `OMP_NUM_THREADS=1` is set, Eigen's OpenMP-backed parallelism is already constrained; the explicit `setNbThreads(1)` is belt-and-braces for any non-OMP Eigen threading paths.

**xdist oversubscription guard:** `benchmarks/conftest.py` hard-asserts `pytest-xdist` is inactive (e.g. `assert "PYTEST_XDIST_WORKER" not in os.environ`). Running benchmarks under `-n auto` multiplies the OMP thread budget across workers, producing noisy and non-reproducible numbers. Fail fast rather than measure garbage.

Baseline JSON is refreshed explicitly (PR-gated — see workflow above), not auto-updated on green builds. This keeps the baseline meaningful as a regression anchor rather than drifting with noise.

### Component: KD-tree parallelization (`9f5`)

`NeighbourFinder<VecMatType>::update()` becomes `#pragma omp parallel for` over the query loop at `NeighbourFinder.hpp:157`. Requires:

- Making `queriedFeature`, `neighbourIndices`, `neighbourSquaredDistances`, and the `knnResultSet` per-thread-local — **declared inside the loop body**, not outside. Current code has `knnResultSet` declared once at `NeighbourFinder.hpp:154` and calls `.init(...)` at top of each iteration (line 159); under parallel-for, that shared result set races on both `init()` and the nanoflann `addPoint` callback. The fix pattern:
  ```cpp
  // WRONG (current pattern, race under parallel for):
  nanoflann::KNNResultSet<FeatureType> knnResultSet(_numNeighbours);
  for (size_t i = 0; i < N; ++i) {
      knnResultSet.init(indices, dists);           // races
      _kdTree->findNeighbors(knnResultSet, ...);
  }

  // RIGHT (declare-inside-loop, thread-safe):
  #pragma omp parallel for
  for (size_t i = 0; i < N; ++i) {
      FeatureType queriedFeature = ...;
      std::vector<size_t> neighbourIndices(_numNeighbours);
      std::vector<FeatureType> neighbourSquaredDistances(_numNeighbours);
      nanoflann::KNNResultSet<FeatureType> knnResultSet(_numNeighbours);
      knnResultSet.init(neighbourIndices.data(), neighbourSquaredDistances.data());
      _kdTree->findNeighbors(knnResultSet, ...);
      // write output at row i (distinct per thread)
  }
  ```
- Verifying `nanoflann`'s `findNeighbors` is thread-safe on a shared `_kdTree`. Confirmed at `vendor/nanoflann.hpp:1243-1256`: `findNeighbors` is `const` with a stack-local `distance_vector_t dists`; mutation happens on stack-local state. `KDTreeSingleIndexAdaptor` is read-only-safe **iff the index is not being rebuilt**. The kd-tree build in `set_source_points` (invoked by `ViscoElasticTransformer::_update_neighbours` on each outer ICP iteration) stays sequential and happens OUTSIDE the parallel region — only the per-query loop is parallelized. `SearchParams(32, 0.0001, true)` with `sorted=true` is safe under parallel queries (sorting is per-query result set, not global).
- Output arrays `_outNeighbourIndices` and `_outNeighbourSquaredDistances` are written at distinct row indices `i` per thread — no race.

**CMake + OpenMP wiring (MSVC-specific):** MSVC's default `/openmp` flag is OpenMP 2.0 and rejects `size_t` loop counters. `library/CMakeLists.txt` must force `/openmp:llvm` AND avoid dragging in the legacy OpenMP 2.0 runtime via `OpenMP::OpenMP_CXX` (which `find_package(OpenMP)` on MSVC resolves to `vcomp`, not the LLVM OpenMP used by `/openmp:llvm`). Committed approach: use `target_compile_options(meshmonk_lib PRIVATE /openmp:llvm)` on MSVC and skip linking `OpenMP::OpenMP_CXX` on that platform (the `/openmp:llvm` flag brings in `libomp` implicitly). On GCC/Clang/Apple Clang, the standard `target_link_libraries(meshmonk_lib PUBLIC OpenMP::OpenMP_CXX)` from ADR D2 applies unchanged:

```cmake
find_package(OpenMP REQUIRED)
if(MSVC)
    target_compile_options(meshmonk_lib PRIVATE /openmp:llvm)
else()
    target_link_libraries(meshmonk_lib PUBLIC OpenMP::OpenMP_CXX)
endif()
```

This is a required default, not a fallback. OpenMesh (vendored subdirectory) does not enable OpenMP on its own targets, so there's no risk of dual-runtime linkage from that side.

Expected: 4–8× speedup on correspondence phase on 8-core machines. Verify against baseline before declaring done.

### Component: Cache temporaries (`6a5`)

Hoist repeated allocations out of the **inner ICP loops** and reuse via `setZero()` / `conservativeResize()`:

- `NonrigidRegistration.cpp:59-61` — `correspondingFeatures`, etc. currently allocated per iteration of `NonrigidRegistration`'s own loop.
- `ViscoElasticTransformer.cpp` — `regularizedForceField`, `unregulatedDisplacementField`, `temporaryDisplacementField` (lines 162, 205, 248), allocated per `update()` call inside `ViscoElasticTransformer`.
- `InlierDetector.cpp:86` — `tempInlierWeights`, allocated per `update()` call inside `InlierDetector`.

Move these to member fields on the respective classes, sized in `set_input`/`set_output` (matching the existing pattern in `ViscoElasticTransformer::set_output:25-27`, which already resizes `_displacementField`/`_oldDisplacementField` on each set), and zeroed at top of each iteration. **Do NOT** size-on-first-call: that pattern breaks if a registration instance is reused with a different-sized floating mesh. Size-in-setter is the robust pattern. A small unit test exercising back-to-back `update()` calls with different input sizes catches regressions.

**Scope boundary — what `6a5` does NOT change.** `NonrigidRegistration::update()` at `library/src/NonrigidRegistration.cpp:83,91` constructs `InlierDetector` and `ViscoElasticTransformer` as **local stack objects**. `PyramidNonrigidRegistration::update()` at `PyramidNonrigidRegistration.cpp:174` constructs a fresh `NonrigidRegistration` per pyramid layer. This means the member-field buffers we hoist are destroyed and re-allocated on every outer call — savings accumulate only within a single `update()` invocation's inner ICP loop, not across pyramid layers or across standalone nonrigid calls. Promoting the owning objects to members of the caller (so they persist across calls) would be a deeper structural change that touches constructor/destructor lifetimes and pyramid layer semantics; per ADR D7 (no scope creep), we keep `6a5` to within-call hoisting and accept the smaller win. If profiling after `6a5` merges shows that outer-call allocations dominate, a follow-up bead can promote the owning classes.

Expected: **3–8% wall-clock on pyramid** (revised down from an initial 5–15% estimate after realizing outer-loop re-construction caps the savings). The inner-ICP allocation savings are real but bounded. `6a5`'s D4 merge threshold is set accordingly (see ADR).

No API change. Numeric output should be bitwise identical (same operations, different memory) — this is enforceable as a merge gate and ADR D4 requires it for `6a5` specifically.

### Component: Smoothing parallelization (`bdt`)

Three loops, all with independent per-vertex work (`#pragma omp parallel for schedule(static)` with chunk size ≥ 16 rows to mitigate false-sharing on contiguous `Vec3Mat` member buffers — each row is 12 bytes, 16 rows ≈ 192 bytes clears a 64-byte cache line comfortably):
- `ViscoElasticTransformer::_update_viscously` smoothing loop (`:167`)
- `ViscoElasticTransformer::_update_elastically` inner loop (`:214`)
- `InlierDetector::_smooth_inlier_weights` inner loop (`:95`)

Each loop reads from one buffer and writes to a different buffer at row `i` — no aliasing. **Invariant (code-review-enforced):** the read-source buffer and the write-destination buffer under 6a5+bdt remain distinct matrices. The existing code is already careful about this (`forceField` → `regularizedForceField`, `_displacementField` → `unregulatedDisplacementField`); after hoisting to member fields in 6a5, both endpoints remain distinct — maintain that. Thread-local scalars (`vectorAverage`, `sumWeights`) stay loop-local by declaration inside the body.

The outlier-diffusion loop (`ViscoElasticTransformer::_update_outlier_transformation`, `:258`) is structurally parallelizable but **excluded from `bdt`'s initial scope** — it runs fewer iterations (`_outlierDiffusionIterations`) and has an inlier-skip branch that reduces per-row work, so thread-launch overhead may exceed the speedup. Per ADR D1 (measure before optimize), this exclusion is not an a-priori decision: `bdt` must take one 30-second measurement against the baseline harness and record the outcome in the bead notes — either include the loop if it shows a ≥1.3× win at 5K, or exclude it with the number that justifies the exclusion. "Likely not worth it" without a measurement is the anti-pattern D1 forbids.

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
- `50K`-vertex scenario uses `data/Template.obj` (if it has ≥50K vertices — verified in `44m`, adjusted if not). We don't have a >100K-vertex benchmark yet; revisit if large-mesh perf becomes a user complaint.
- OpenMP on Apple Silicon macOS with Accelerate-linked Eigen can interact oddly with Grand Central Dispatch. Benchmark on macOS before declaring `9f5` done.
- Calling meshmonk register functions from multiple Python threads concurrently (e.g. a `ThreadPoolExecutor` batching registrations) oversubscribes OpenMP — each Python worker spawns `N_CPU` OMP threads. Not a correctness bug, but a perf footgun. Users running batch workloads should set `OMP_NUM_THREADS=1` and parallelize at the Python level instead. Documented in `benchmarks/README.md` and `docs/gotchas.md` (once `gw1` lands).
- `benchmarks/baseline.json` must be refreshed on GitHub Actions runner-image upgrades (`ubuntu-latest` rolls roughly yearly) OR every minor release, whichever comes first. Stale baselines turn the regression gate into noise. Track via a scheduled reminder bead.

---

## Open Questions

1. **Baseline refresh cadence** — per minor release, per runner-image upgrade, or both? Known Limitations commits to "whichever comes first"; revisit after a week of CI noise data from 44m.
2. **Should `6a5` land before or after `9f5`?** Both depend only on `44m`. `6a5` is riskier in terms of subtle state bugs (hoisted buffers might leak state across calls on re-entry) but lower in raw speedup. Recommend `9f5` first — bigger win, simpler to review.

