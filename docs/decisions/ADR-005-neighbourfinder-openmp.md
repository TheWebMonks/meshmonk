# ADR-005: NeighbourFinder OpenMP Parallelization

**Status:** Accepted
**Date:** 2026-04-23
**Design:** [NeighbourFinder OpenMP design doc](../../history/2026-04-23-neighbourfinder-openmp-design.md)
**Parent:** [ADR-001 meshmonk-modernization](./ADR-001-meshmonk-modernization.md)
**Evidence:** [Hotspot profile 20260423](../../docs/perf/hotspot-profile-20260423-0e3072f.md)

## Context

The hotspot profile (748 profiling bead) identified the correspondence pipeline
as the dominant cost: `SymmetricCorrespondenceFilter` (22%), `CorrespondenceFilter`
(16%), and `NeighbourFinder` (15%) together account for ~53% of nonrigid
wall-time at 100K vertices with super-linear scaling (k≈1.23–1.26). These three
labels are a composition stack — `NeighbourFinder` is the leaf; parallelizing it
fixes all three simultaneously.

The immediate cause is a sequential for-loop over `_numQueriedElements` vertices
in `NeighbourFinder::update()`, each performing an independent k-NN query
against a shared read-only KD-tree.

## Decisions

### D1: Parallelize at NeighbourFinder, not at CorrespondenceFilter or SymmetricCorrespondenceFilter level

**Firmness: FIRM**

The `#pragma omp parallel for` goes on the inner query loop in `NeighbourFinder::update()`,
not on the push/pull update calls in `SymmetricCorrespondenceFilter`.

**Rationale:** A single change at `NeighbourFinder` propagates through the full
stack: `CorrespondenceFilter` (which holds one `NeighbourFinder`) and
`SymmetricCorrespondenceFilter` (which wraps two `CorrespondenceFilter`s) both
benefit automatically. Parallelizing at `SymmetricCorrespondenceFilter` level
(running push and pull as concurrent tasks) would require OpenMP tasks or
sections, only halves the symmetric cost, and does not help the
`CorrespondenceFilter` non-symmetric path. The inner-loop site is embarrassingly
parallel; the task-level site is not.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **NeighbourFinder inner loop (chosen)** | Single change site, fixes all three labels, embarrassingly parallel | — |
| SymmetricCorrespondenceFilter push+pull as tasks | No change to NeighbourFinder | Only 2× on symmetric path; doesn't help CorrespondenceFilter; more complex (task synchronization) |
| Affinity construction loop | Also parallelizable | ~10% of correspondence cost — not the bottleneck |

**What would invalidate this:** If the profiling characterization turns out to be
wrong — e.g., affinity construction dominates at a scale tier not yet measured —
adding parallelization there too would be warranted.

---

### D2: Thread-local buffers via in-loop declaration (not `thread_local` or pre-allocated pools)

**Firmness: FIRM**

The local buffers (`queriedFeature`, `neighbourIndices`, `neighbourSquaredDistances`,
`knnResultSet`) are moved from before the loop to inside the loop body. Each
OpenMP thread allocates them on its stack per iteration.

**Rationale:** In-loop declaration is the minimal correct change: it makes the
buffers automatically per-thread without any machinery. The allocation cost
(small `std::vector` per iteration) is negligible relative to the k-NN query
time (ms range at 10K+). `thread_local` would give persistent storage across
calls but `KNNResultSet` is stateful per-query (requires `init()` before each
use), so persistent storage would need explicit reinitialization anyway — no win.
Pre-allocated pools add complexity with no measurable benefit.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **In-loop declaration (chosen)** | Minimal change, zero machinery, correct by construction | Per-iteration allocation (negligible at ms-scale per-call) |
| `thread_local` static buffers | Avoids per-iteration allocation | `KNNResultSet` must be reinit'd each query anyway; adds lifecycle complexity |
| Pre-allocated per-thread pools | No allocation in hot path | Significant added machinery; premature optimization |

**What would invalidate this:** Profiling showing that vector allocation overhead
is measurable (would appear as significant overhead at 1K tier where per-call
time is ~0.2ms). If so, switch to `thread_local` with explicit `resize` +
`init`.

---

### D3: `schedule(static)` for the parallel for

**Firmness: FLEXIBLE**

The OpenMP schedule is `static` (default chunk: n_vertices / n_threads).

**Rationale:** Meshes in the production workload have roughly uniform vertex
density, so k-NN query cost per vertex is approximately constant. Static
scheduling has lower overhead than `dynamic` (no runtime work-stealing) and
achieves near-equal load balance for uniform work. `dynamic` adds per-chunk
queue overhead, which matters when per-iteration cost is small (1K tier,
~0.12µs/query).

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **`schedule(static)` (chosen)** | Zero runtime overhead, near-equal balance for uniform meshes | Sub-optimal if vertex density is highly non-uniform |
| `schedule(dynamic, 64)` | Better load balance for non-uniform meshes | Per-chunk queue overhead; measurable at 1K tier |
| `schedule(guided)` | Adaptive chunk size | More complex, rarely better than static for k-NN |

**What would invalidate this:** A use case with high local-density variation (e.g.
a mesh with clustered high-resolution regions). Switch to `schedule(dynamic, 64)`
and verify with the profiling driver.

---

### D4: Guarded `find_package(OpenMP)` — silent no-op when absent

**Firmness: FIRM**

`library/CMakeLists.txt` adds:

```cmake
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(meshmonk_lib PUBLIC OpenMP::OpenMP_CXX)
endif()
```

When OpenMP is not found (e.g. macOS Apple Clang without libomp), the build
succeeds and `NeighbourFinder::update()` runs serially — the `#pragma omp`
directive is silently ignored by compilers that don't recognize it (C++ compilers
treat unknown pragmas as warnings, not errors).

**Rationale:** OpenMP is not available by default on macOS Apple Clang. An
unconditional `find_package(OpenMP REQUIRED)` would break developer builds on
macOS out of the box. The guard keeps the library buildable everywhere; CI on
Linux (GCC) links OpenMP and measures the speedup; macOS developers install
libomp via Homebrew when they want parallel performance.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **Guarded find_package (chosen)** | Builds everywhere, parallel on GCC/Linux out of the box | Serial on macOS without libomp |
| `find_package(OpenMP REQUIRED)` | Guarantees parallel build | Breaks macOS Apple Clang builds by default |
| `MESHMONK_OPENMP` build option | Explicit opt-in | Parallel should be the default when available; opt-in creates a permanently-serial default |

**What would invalidate this:** Requiring parallel builds on macOS as a CI
target — at that point the Homebrew libomp path should be documented and added
to the macOS CI workflow.

---

### D5: KD-tree build remains serial

**Firmness: FIRM**

`set_source_points()` (which calls `buildIndex()`) is not parallelized. It is
called once per ICP iteration before the query loop and is not in the measured
hotspot.

**Rationale:** The profiling data shows no measurable cost for tree construction
— it does not appear as a significant label. Parallelizing it would require
either a parallel KD-tree library or splitting the point set across threads, both
of which are significant complexity increases with zero demonstrated ROI.

**What would invalidate this:** Tree build time becoming visible in the hotspot
profile at a new scale tier (e.g. 500K+ vertices). At that point a parallel
build via a dedicated library (e.g. nanoflann's own parallel build support, if
added upstream) is the right path.
