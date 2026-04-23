# Profiling Infrastructure Design

**Status:** Accepted  
**Date:** 2026-04-23  
**Decisions:** [ADR-004](../docs/decisions/ADR-004-profiling.md)  
**Bead:** `meshmonk-modernization-748`

## What We're Building

Scale-aware profiling infrastructure to map where MeshMonk registration actually spends its time at 1K / 10K / 100K vertex meshes. The output is a committed hotspot report that gates every future perf bead — no optimization is implemented here.

## Mesh Acquisition

Template.obj (7160 verts) and demoFace.obj are tessellated via loop subdivision (trimesh) and then downsampled to hit the target scale:

- **1K**: downsample both at runtime via existing `_harness.downsample_to`
- **10K**: one loop subdivision (~28K) → downsample to 10K; commit as `data/Template_10K.obj` + `data/DemoFace_10K.obj`
- **100K**: two loop subdivisions (~112K) → downsample to 100K; commit as `data/Template_100K.obj` + `data/DemoFace_100K.obj`

Generation script: `scripts/generate_profiling_meshes.py` — run once, not part of CI.

Both meshes matched at each scale tier so k-NN scaling curves have one variable (n). Template+DemoFace pairs preserve realistic shape difference vs. registering a mesh against itself.

## RAII Timer Infrastructure

A single header `library/include/meshmonk/profiling.hpp`, active only when `-DMESHMONK_PROFILING=ON`.

- `ScopedTimer` RAII guard: records elapsed microseconds into a global `g_profiler` accumulator on destruction
- `g_profiler` accumulates `{label → {total_us, count}}` across a full registration run
- Callsite pattern: `#ifdef MESHMONK_PROFILING auto _t = g_profiler.scoped("NeighbourFinder::update"); #endif`

Labels to instrument (initial set — expand during implementation if callsites reveal additional cost):

| Label | Location |
|---|---|
| `NeighbourFinder::update` | NeighbourFinder.cpp |
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

## Build Integration

`library/CMakeLists.txt`:

```cmake
option(MESHMONK_PROFILING "Enable registration profiling instrumentation" OFF)
if(MESHMONK_PROFILING)
  target_compile_definitions(meshmonk_lib PUBLIC MESHMONK_PROFILING)
endif()
```

Same option propagated to `_meshmonk_core` in the root `CMakeLists.txt`. Default OFF — CI and baseline benchmarks are unaffected.

## Python Bindings Surface

Three functions always compiled into the extension (no-ops when `MESHMONK_PROFILING` not set):

- `meshmonk.profiling_reset()` — clear accumulator
- `meshmonk.profiling_dump() → dict` — return `{label: {total_us, count}}`, reset accumulator
- `meshmonk.profiling_enabled() → bool` — check build flag

## Profiling Driver

`profiling/run_profile.py` — standalone script, not part of CI. For each scale tier (1K / 10K / 100K) and each registration mode (rigid / nonrigid / pyramid):

1. `profiling_reset()`
2. Run registration
3. `data = profiling_dump()`
4. Compute per-label: `ms_per_call = total_us / count / 1000`, `share = total_us / run_total_us`, `amdahl_ceiling = 1 / (1 - share)`
5. Fit scaling exponent k per label via log-log regression across the three tiers

## Report Format

`docs/perf/hotspot-profile-YYYYMMDD.md` — committed after each profiling run.

Sections:
1. **Methodology** — mesh provenance, build config, machine spec, run command
2. **Per-tier ranked table** — columns: label, share (%), ms/call, invocations, Amdahl ceiling
3. **Scaling exponents** — per-label k from log-log fit, flagged as O(n), O(n log n), O(n²)
4. **Recommendation list** — paths worth a future perf bead: `share > 5%` AND `ms/call > 0.03ms` (OpenMP-viable regime); each recommendation includes suggested tool (OpenMP / algorithmic / Eigen-level) and the arithmetic for why it can clear the ≥5% gate

Paths below the granularity threshold are listed with a note explaining why OpenMP is not viable (fork-join overhead lesson from 9f5/bdt).
