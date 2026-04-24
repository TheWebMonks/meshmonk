# NeighbourFinder Line-Level Profiling Design

**Status:** Accepted
**Date:** 2026-04-23
**Bead:** `meshmonk-modernization-5jf`
**Parent infra:** [ADR-004 profiling](../docs/decisions/ADR-004-profiling.md)
**Decisions:** [ADR-006](../docs/decisions/ADR-006-neighbourfinder-profiling.md)

## What We're Building

A line-level cost breakdown of `NeighbourFinder::update()` at 100K vertices using the existing RAII `ScopedTimer` infrastructure. The output is a committed profiling report that identifies which of four cost buckets dominates, gating the choice of the next optimization bead.

## Context

ADR-005 (OpenMP parallelization) was closed without implementation — prior bead `9f5` measured a regression (0.93–0.99×) at 7K with the identical `#pragma omp parallel for`, and the Amdahl ceiling is ~1.17× at 100K (~15% leaf share). Before committing to any optimization direction, we need to know WHERE inside `NeighbourFinder::update()` the cost lives.

The four cost buckets inside the query loop:

| Bucket | Code | Intervention if dominant |
|---|---|---|
| `buffer_alloc` | `std::vector` construction for `queriedFeature`, `neighbourIndices`, `neighbourSquaredDistances` | `thread_local` storage (revisit ADR-005 D2) |
| `query_setup` | `KNNResultSet.init()` + feature extraction from Eigen row | Per-query setup reduction |
| `tree_query` | `_kdTree->index->findNeighbors(...)` | Neighbour caching (#1) or SIMD batching (#2) |
| `result_copy` | Copying indices/distances to output matrices | Vectorize or eliminate copy |

## Instrumentation Approach

Subdivide the query loop in `NeighbourFinder::update()` using the existing `ScopedTimer` RAII pattern. Each `ScopedTimer` scope wraps exactly one cost bucket. Labels follow the existing naming convention from ADR-004 D2.

```cpp
for (int i = 0; i < _numQueriedElements; i++) {
    #ifdef MESHMONK_PROFILING
    auto _t_alloc = g_profiler.scoped("NeighbourFinder::buffer_alloc");
    #endif
    std::vector<float> queriedFeature(_numDimensions);
    std::vector<size_t> neighbourIndices(_numNeighbours);
    std::vector<float> neighbourSquaredDistances(_numNeighbours);
    nanoflann::KNNResultSet<float> knnResultSet(_numNeighbours);
    #ifdef MESHMONK_PROFILING
    _t_alloc.~ScopedTimer();  // explicit early end — NOT the right pattern; see D1
    #endif
    ...
}
```

**D1 (FIRM): Use RAII scope bracketing, not explicit destructor calls.** Introduce a nested scope `{}` for each bucket to trigger RAII destruction at the right point. This is correct C++ and avoids UB from explicit destructor calls.

Correct pattern:
```cpp
for (int i = 0; i < _numQueriedElements; i++) {
    {
        #ifdef MESHMONK_PROFILING
        auto _t = g_profiler.scoped("NeighbourFinder::buffer_alloc");
        #endif
        queriedFeature.resize(_numDimensions);
        neighbourIndices.resize(_numNeighbours);
        neighbourSquaredDistances.resize(_numNeighbours);
        knnResultSet = nanoflann::KNNResultSet<float>(_numNeighbours);
    }
    {
        #ifdef MESHMONK_PROFILING
        auto _t = g_profiler.scoped("NeighbourFinder::query_setup");
        #endif
        // feature extraction + knnResultSet.init(...)
    }
    {
        #ifdef MESHMONK_PROFILING
        auto _t = g_profiler.scoped("NeighbourFinder::tree_query");
        #endif
        _kdTree->index->findNeighbors(...);
    }
    {
        #ifdef MESHMONK_PROFILING
        auto _t = g_profiler.scoped("NeighbourFinder::result_copy");
        #endif
        // copy indices/distances to output matrices
    }
}
```

This requires hoisting buffer declarations above the loop. Timers must remain outside `#pragma omp parallel` regions (ADR-004 D2 thread-safety note — no OpenMP here so this is not a concern).

**D2 (FIRM): Hoist loop buffers above the loop.** Move `queriedFeature`, `neighbourIndices`, `neighbourSquaredDistances`, `knnResultSet` declarations before the `for` loop, `resize()`/`reinit()` inside the loop. This separates first-allocation cost (amortized) from per-iteration cost (what we're measuring).

**D3 (FIRM): Run at 100K nonrigid only.** The existing per-tier infrastructure is not needed for this investigation — we only need to identify which bucket dominates, which is visible at the tier where NeighbourFinder is measurably hot. The driver runs one nonrigid registration at 100K with ≥5 measured repetitions (1 warmup discarded), median across reps.

## Files to Change

- `library/include/meshmonk/NeighbourFinder.hpp` — add 4 ScopedTimer scopes inside the query loop; hoist buffer declarations above the loop

## Build

```bash
pip install -C cmake.args=-DMESHMONK_PROFILING=ON -e .
```

## Driver Invocation

```python
python -m profiling.run_profile \
  --tiers 100k \
  --modes nonrigid \
  --runs 5 \
  --warmup 1 \
  --seed 42 \
  --out docs/perf/neighbourfinder-breakdown-20260423-<shortsha>.md
```

## Output Artifact

`docs/perf/neighbourfinder-breakdown-YYYYMMDD-<shortsha>.md` committed to the repo. Sections:
1. **Methodology** — host, compiler, build flags, git SHA, run command
2. **Per-bucket table** — label, share (%), ms/call (overhead-corrected), invocations
3. **Recommendation** — which next bead to create based on the dominant bucket (see bead 5jf notes for the decision tree)
