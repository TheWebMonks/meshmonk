# Profiling Infrastructure Design

**Status:** Accepted
**Date:** 2026-04-23
**Decisions:** [ADR-004](../docs/decisions/ADR-004-profiling.md)
**Bead:** `meshmonk-modernization-748`

## What We're Building

Scale-aware profiling infrastructure to map where MeshMonk registration actually spends its time at 1K / 10K / 100K vertex meshes. The output is a committed hotspot report that gates every future perf bead — no optimization is implemented here.

## Mesh Acquisition

All three tiers share one subdivision-then-downsample pedigree so k-NN cost has exactly one variable per tier (n). Template.obj (7160 verts) and demoFace.obj are the seeds.

- **1K**: one loop subdivision (~28K) → downsample to 1K; commit as `data/Template_1K.obj` + `data/DemoFace_1K.obj`
- **10K**: one loop subdivision (~28K) → downsample to 10K; commit as `data/Template_10K.obj` + `data/DemoFace_10K.obj`
- **100K**: two loop subdivisions (~112K) → downsample to 100K; commit as `data/Template_100K.obj` + `data/DemoFace_100K.obj`

Generation script: `scripts/generate_profiling_meshes.py` — run once, not part of CI. Mesh file SHAs are recorded in the report Methodology so a regeneration invalidates old reports.

Both meshes matched at each scale tier so k-NN scaling curves have one variable (n). Template+DemoFace pairs preserve realistic shape difference vs. registering a mesh against itself.

## RAII Timer Infrastructure

A single header `library/include/meshmonk/profiling.hpp`, active only when `-DMESHMONK_PROFILING=ON`.

- `ScopedTimer` RAII guard: records elapsed microseconds into a global `g_profiler` accumulator on destruction
- `g_profiler` is declared `inline` in the header (C++20 baseline per ADR-001 D1), so a single definition survives across all translation units with no `.cpp` needed
- `g_profiler` accumulates `{label → {total_us, count}}` across a full registration run
- Callsite pattern: `#ifdef MESHMONK_PROFILING auto _t = g_profiler.scoped("NeighbourFinder::update"); #endif`

**Thread safety:** labels must be placed outside `#pragma omp parallel` regions. The accumulator is single-threaded today (matching ADR-002 D1's `OMP_NUM_THREADS=1` baseline convention, which the profiling driver also pins). If per-thread profiling is needed later, switch the accumulator to `thread_local` storage with a merge step on dump — out of scope for v0.1.

**Timer overhead calibration:** implementation includes a calibration step that measures empty-`ScopedTimer` overhead on the target host (expected 60–160ns per scope from two `std::chrono::steady_clock::now()` calls). The per-host overhead is recorded in the report Methodology and the driver subtracts `count × overhead` from each label's `total_us` before computing share and `ms/call`. Any label with `ms/call < 5 × timer_overhead` is marked **UNRELIABLE** and excluded from the recommendation list (this matters at the 30µs OpenMP-viability threshold, where uncorrected timer overhead alone is ~0.3–0.5% of the measured value).

**Pyramid-mode label scoping:** `PyramidNonrigidRegistration::update()` invokes `NonrigidRegistration::update()` once per pyramid layer at different downsampled mesh sizes. To keep per-label stats and scaling fits meaningful, pyramid callsites append a `/layer{N}` suffix (e.g. `NonrigidRegistration::update/layer0`, `NonrigidRegistration::update/layer1`). The standalone nonrigid mode uses the unsuffixed label. The report dedupes across modes in the ranking and fit.

Labels to instrument (initial set — expand during implementation if callsites reveal additional cost):

| Label | Location |
|---|---|
| `NeighbourFinder::update` | NeighbourFinder.hpp (template method body, ~line 144) |
| `CorrespondenceFilter::update` | CorrespondenceFilter.cpp |
| `SymmetricCorrespondenceFilter::update` | SymmetricCorrespondenceFilter.cpp |
| `InlierDetector::update` | InlierDetector.cpp |
| `ViscoElasticTransformer::_update_viscously` | ViscoElasticTransformer.cpp |
| `ViscoElasticTransformer::_update_elastically` | ViscoElasticTransformer.cpp |
| `ViscoElasticTransformer::_update_outlier_transformation` | ViscoElasticTransformer.cpp |
| `ViscoElasticTransformer::_update_smoothing_weights` | ViscoElasticTransformer.cpp |
| `ViscoElasticTransformer::_apply_transformation` | ViscoElasticTransformer.cpp |
| `NonrigidRegistration::update` (full ICP loop) | NonrigidRegistration.cpp |
| `RigidRegistration::update` | RigidRegistration.cpp |
| `PyramidNonrigidRegistration::update` | PyramidNonrigidRegistration.cpp |

Pyramid-mode callsites append `/layer{N}` to the `NonrigidRegistration::update` label as described above.

## Build Integration

`library/CMakeLists.txt`:

```cmake
option(MESHMONK_PROFILING "Enable registration profiling instrumentation" OFF)
if(MESHMONK_PROFILING)
  target_compile_definitions(meshmonk_lib PUBLIC MESHMONK_PROFILING)
endif()
```

The `PUBLIC` scope on `meshmonk_lib` propagates `MESHMONK_PROFILING` transitively to `_meshmonk_core` via its `PRIVATE` link (root `CMakeLists.txt:46`) — no additional `target_compile_definitions` on `_meshmonk_core` is needed. Root `CMakeLists.txt` declares the same `option(...)` so `pip install -C cmake.args=-DMESHMONK_PROFILING=ON .` (scikit-build-core pass-through per ADR-001 D5) works for a profiling build. Default OFF — CI and baseline benchmarks are unaffected.

## Python Bindings Surface

Five functions always compiled into the extension (no-ops when `MESHMONK_PROFILING` not set):

- `meshmonk.profiling_reset()` — clear accumulator
- `meshmonk.profiling_peek() → dict` — return `{label: {total_us, count}}` without reset; safe to call multiple times
- `meshmonk.profiling_dump() → dict` — return `{label: {total_us, count}}` AND reset accumulator (equivalent to `peek()` then `reset()`); the reset side-effect is in the name for clarity
- `meshmonk.profiling_enabled() → bool` — check build flag
- `meshmonk.profiling_calibrate(n: int = 1_000_000) → int` — run n empty `ScopedTimer` scopes and return measured nanoseconds per scope (returns 0 when built without `MESHMONK_PROFILING`)

The `peek` / `dump` split prevents a common footgun: if a mid-registration progress print calls `dump()`, the final `dump()` sees near-empty data. The driver uses `reset → run → dump` (single-dump contract). Ad-hoc callers use `peek` for read-only snapshots. The `calibrate` function measures true empty-`ScopedTimer` overhead in C++ so the driver doesn't have to estimate it from Python-side timings.

## Profiling Driver

`profiling/run_profile.py` — standalone script, not part of CI.

**Step 0:** `assert meshmonk.profiling_enabled()`; if False, print clear error ("meshmonk was not built with -DMESHMONK_PROFILING=ON — the report would be all zeros") and exit 2. This prevents a silent failure mode that would otherwise surface only when a reader notices the zeros.

**Step 1 (calibration):** measure empty-`ScopedTimer` overhead via a loop of N=10⁶ no-op scopes; record per-scope nanoseconds in the report Methodology and use it to correct measured totals.

**For each scale tier (1K / 10K / 100K) and each registration mode (rigid / nonrigid / pyramid):**

1. `profiling_reset()`
2. Run registration (one run per repetition)
3. `data = profiling_dump()`
4. Repeat steps 1–3 **≥5 times** per (tier, mode); one warmup run is discarded before the 5 measured runs
5. For each label, take the median `total_us` and median `count` across the 5 runs; compute per-label: `ms_per_call = (total_us - count × overhead_us) / count / 1000`, `share = total_us_corrected / run_total_us_corrected`, `amdahl_ceiling = 1 / (1 - share)`, `expected_fork_join_cost_us = 15 × count` (for the OpenMP recommendation logic)
6. Fit scaling exponent `k` per label via log-log regression across the three tiers. **Because three points yield zero residual degrees of freedom, the fit is reported as indicative only.** A fourth tier (e.g. 30K) or confirmation via per-tier variance is required before any perf bead cites the exponent as load-bearing.

### Running

```
python -m profiling.run_profile \
  --tiers 1k,10k,100k \
  --modes rigid,nonrigid,pyramid \
  --runs 5 \
  --warmup 1 \
  --seed 42 \
  --out docs/perf/hotspot-profile-YYYYMMDD-<shortsha>.md
```

`OMP_NUM_THREADS=1` and `EIGEN_DONT_PARALLELIZE=1` are pinned by the driver to match ADR-002 D1 baseline convention. Seed 42 matches the benchmark harness convention (`bench_pyramid.py`).

## Report Format

`docs/perf/hotspot-profile-YYYYMMDD-<shortsha>.md` — committed after each profiling run. The short git SHA disambiguates multiple runs per day and encodes build provenance.

Sections:

1. **Methodology** — enumerated fields:
   - CPU model + core count
   - OS + kernel version
   - Compiler + version
   - `CMAKE_BUILD_TYPE`, `-O` level, LTO on/off
   - OpenMesh / nanoflann / Eigen versions (vendored — pinned by git SHA)
   - Repo git SHA
   - `OMP_NUM_THREADS`, `EIGEN_DONT_PARALLELIZE` at runtime
   - Mesh file SHA-256 for each tier's Template + DemoFace
   - Run command + seed
   - Empty-`ScopedTimer` overhead (calibrated, ns/scope)
   - Repetitions per (tier, mode) and median-selection policy
2. **Per-tier ranked table** — columns: label, share (%), ms/call (corrected), invocations, Amdahl ceiling, expected fork-join cost (15µs × count), UNRELIABLE flag if ms/call < 5× overhead
3. **Scaling exponents** — per-label `k` from log-log fit across 1K/10K/100K. Reported as indicative exponent with buckets k≈1.0 / k≈1.16 / k≈2.0 as approximate complexity classes; explicitly noted that three-point fits cannot reliably discriminate these buckets and require a fourth tier or per-tier variance before load-bearing use.
4. **Recommendation list** — paths worth a future perf bead: `share > 5%` AND `ms/call > 0.03ms` (OpenMP-viable regime) AND NOT UNRELIABLE; each recommendation includes suggested tool (OpenMP / algorithmic / Eigen-level), Amdahl ceiling arithmetic, and expected fork-join cost (so a small-share path with high invocation count is visibly disqualified from OpenMP before any bead cites it).

Paths below the granularity threshold are listed with a note explaining why OpenMP is not viable (fork-join overhead lesson from 9f5/bdt). Paths flagged UNRELIABLE are listed in a separate appendix with their corrected ms/call to justify the exclusion.
