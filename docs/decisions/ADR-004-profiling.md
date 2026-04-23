# ADR-004: Profiling Infrastructure — Mesh Acquisition, Timer Design, Build Integration

**Status:** Accepted  
**Date:** 2026-04-23  
**Design:** [Profiling design doc](../../history/2026-04-23-profiling-design.md)  
**Parent:** [ADR-001 meshmonk-modernization](./ADR-001-meshmonk-modernization.md)  
**Related:** Bead `meshmonk-modernization-748`

## Context

Three failed perf beads (6a5, 9f5, bdt) shared a common failure mode: pre-triage estimated a code path looked expensive, implementation proceeded, measured speedup failed the ≥5% gate. ~92% of nonrigid wall-clock at 7K was never measured. This bead builds systematic profiling infrastructure to produce a hotspot map before any further optimization work begins.

## Decisions

### D1: Subdivide + downsample for all tiers (1K / 10K / 100K)

**Firmness: FIRM**

Loop subdivision via trimesh scales Template.obj (7160 verts) and demoFace.obj to ~28K (one pass) and ~112K (two passes), then downsampled to the target tier. **All three tiers share this subdivision-then-downsample pedigree** — the 1K tier is generated via one subdivision pass then downsample to 1K, not by downsampling the raw 7160-vert Template directly. Generated once, committed to `data/` as `Template_1K.obj`, `Template_10K.obj`, `Template_100K.obj`, `DemoFace_1K.obj`, `DemoFace_10K.obj`, `DemoFace_100K.obj`.

**Rationale:** Subdivision preserves mesh topology and produces geometry with realistic face density. No dependency on external scan datasets. Template+DemoFace pairs preserve the shape difference present in production workloads. Both meshes matched at each scale tier so k-NN cost has one variable (n) and scaling exponents are clean. Using a single pedigree for all three tiers (rather than downsampling 1K from the raw Template while generating 10K/100K from subdivisions) keeps the "one variable per tier" claim honest — mixing pedigrees would change local density distribution at 1K and produce a spurious scaling-curve discontinuity.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **Subdivide + downsample (chosen)** | Realistic topology, one variable per tier, no external deps | Subdivided geometry may not match production density exactly |
| Real scan acquisition | Most realistic | External dependency, time-consuming to source and clean |
| Use same mesh for float and target | Simpler (one OBJ per tier) | Removes shape difference — less realistic registration scenario |

**What would invalidate this:** A downstream user reporting that the profiling results don't match their production mesh density — at that point real scan acquisition is warranted.

---

### D2: RAII scoped timer with global accumulator

**Firmness: FIRM**

A single header `library/include/meshmonk/profiling.hpp` defines `ScopedTimer` (RAII guard) and a global `g_profiler` that accumulates `{label → {total_us, count}}` across a full registration run. Active only under `#ifdef MESHMONK_PROFILING`.

**Rationale:** RAII ensures timers are always closed even on early returns. The global accumulator lets the Python driver call `profiling_dump()` once after registration completes rather than managing per-callsite output. The `#ifdef` guard means exactly zero overhead when not profiling — no branch, no dead code surviving the compiler.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **RAII + global accumulator (chosen)** | Zero overhead when off, always closed, single dump point | Global state requires reset between runs |
| Per-callsite file write | No global state | I/O per callsite adds noise to timing; can't aggregate |
| Env-var runtime toggle | Toggle without recompile | Branch-check overhead per callsite; measurable at µs scale |

**What would invalidate this:** Thread-safety requirements — if registration is ever parallelised at the registration-object level, the global accumulator needs a mutex or thread-local storage.

---

### D3: Build flag `-DMESHMONK_PROFILING=ON`, default OFF

**Firmness: FIRM**

CMake `option(MESHMONK_PROFILING ... OFF)` in `library/CMakeLists.txt`. Propagated to `_meshmonk_core` via root `CMakeLists.txt`. CI, `baseline.json`, and the existing benchmark harness are entirely unaffected.

**Rationale:** Compile-time elimination guarantees zero overhead. The recompile cost is irrelevant — profiling always requires a dedicated build (different optimisation level, no inlining noise). An env-var approach would leave branch-check overhead on every instrumented callsite.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **Build flag, default OFF (chosen)** | Guaranteed zero overhead, standard CMake pattern | Requires recompile to toggle |
| Env-var runtime toggle | No recompile | Branch overhead per callsite; contaminates timing |
| Always-on timing | Simpler | Baseline contamination — defeats the purpose |

**What would invalidate this:** A scenario where profiling must be toggled in an installed/distributed build without recompilation — at that point an env-var gate with carefully measured overhead is the fallback.

---

### D4: Five-function Python binding surface with peek/dump split + calibrate

**Firmness: FLEXIBLE**

`profiling_reset()`, `profiling_peek() → dict`, `profiling_dump() → dict`, `profiling_enabled() → bool`, `profiling_calibrate(n: int) → int` are always compiled into the nanobind extension. When `MESHMONK_PROFILING` is not set they are no-ops / return empty / return 0. The Python driver controls the profiling lifecycle.

- `peek()` returns the accumulator snapshot without side-effects (safe to call multiple times)
- `dump()` returns the snapshot AND resets (the reset side-effect is in the name for clarity)
- `calibrate(n)` runs `n` empty `ScopedTimer` scopes in C++ and returns measured ns/scope, so the driver records real empty-timer overhead on the target host rather than a Python-side estimate

**Rationale:** Python-side control keeps the profiling driver simple — reset, run, dump, analyse. No file path coupling between C++ and Python. `profiling_enabled()` lets the driver fail fast with a clear message if the extension was built without the flag. The `peek`/`dump` split prevents a common footgun — a mid-run progress print calling `dump()` would silently consume the accumulator; ad-hoc callers use `peek()` for read-only snapshots while the driver uses the `reset → run → dump` contract.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **peek + dump split (chosen)** | Read-only snapshots are safe; reset semantics named explicitly | Four bindings instead of three |
| Single `dump()` with reset side-effect | Three bindings instead of four | Silent footgun when called twice |
| C++ writes JSON to file directly | No new bindings | File path coupling; harder to test |

**What would invalidate this:** If profiling is needed from a context without Python (e.g. a C++ benchmark binary) — at that point a `dump_to_file()` C++ method is the right addition.
