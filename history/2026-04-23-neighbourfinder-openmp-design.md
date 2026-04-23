# NeighbourFinder OpenMP Parallelization Design

**Status:** Accepted
**Date:** 2026-04-23
**Decisions:** [ADR-005](../docs/decisions/ADR-005-neighbourfinder-openmp.md)
**Parent bead:** `meshmonk-modernization-748` (profiling infrastructure)

## What We're Building

A single-site OpenMP parallelization of the per-vertex k-NN query loop in
`NeighbourFinder::update()`. This is the highest-ROI change identified by the
hotspot profile: one loop change fixes the dominant cost path across all three
correspondence labels simultaneously.

## Why Here

The hotspot report (`docs/perf/hotspot-profile-20260423-0e3072f.md`) shows that
the correspondence pipeline accounts for ~50–65% of nonrigid wall-time at
10K–100K vertices:

| Label | Share (nonrigid, 100K) | k |
|---|---|---|
| `SymmetricCorrespondenceFilter::update` | 22% | 1.23 |
| `CorrespondenceFilter::update` | 16% | 1.23 |
| `NeighbourFinder::update` | 15% | 1.26 |

These three labels are not independent costs — they form a composition stack:

- `NeighbourFinder` — leaf: per-vertex k-NN query loop against a nanoflann KD-tree
- `CorrespondenceFilter` — wraps one `NeighbourFinder` (float→target direction)
- `SymmetricCorrespondenceFilter` — wraps two `CorrespondenceFilter`s (push float→target, pull target→float) and fuses their affinities

`SymmetricCorrespondenceFilter` ≈ 2× `CorrespondenceFilter` ≈ 2× `NeighbourFinder` + affinity
construction overhead (~10% of correspondence cost). Parallelizing `NeighbourFinder::update()`
propagates the speedup through the entire stack.

## The Loop

`library/include/meshmonk/NeighbourFinder.hpp`, line 161:

```cpp
// Current: sequential, shared local buffers
std::vector<float> queriedFeature(_numDimensions);
std::vector<size_t> neighbourIndices(_numNeighbours);
std::vector<float> neighbourSquaredDistances(_numNeighbours);
nanoflann::KNNResultSet<float> knnResultSet(_numNeighbours);

for (; i < _numQueriedElements; ++i) {
    knnResultSet.init(&neighbourIndices[0], &neighbourSquaredDistances[0]);
    for (j = 0; j < _numDimensions; ++j) {
        queriedFeature[j] = (*_inQueriedPoints)(i, j);
    }
    _kdTree->index->findNeighbors(knnResultSet, &queriedFeature[0], ...);
    for (j = 0; j < _numNeighbours; ++j) {
        _outNeighbourIndices(i, j) = neighbourIndices[j];
        _outNeighbourSquaredDistances(i, j) = neighbourSquaredDistances[j];
    }
}
```

Each iteration is an independent k-NN query against a read-only KD-tree. No
data dependency between iterations.

## The Change

Move the local buffers inside the loop body (making them per-thread stack
allocations) and add `#pragma omp parallel for`:

```cpp
#pragma omp parallel for schedule(static)
for (int i = 0; i < (int)_numQueriedElements; ++i) {
    std::vector<float> queriedFeature(_numDimensions);
    std::vector<size_t> neighbourIndices(_numNeighbours);
    std::vector<float> neighbourSquaredDistances(_numNeighbours);
    nanoflann::KNNResultSet<float> knnResultSet(_numNeighbours);

    knnResultSet.init(&neighbourIndices[0], &neighbourSquaredDistances[0]);
    for (int j = 0; j < (int)_numDimensions; ++j) {
        queriedFeature[j] = (*_inQueriedPoints)(i, j);
    }
    _kdTree->index->findNeighbors(knnResultSet, &queriedFeature[0], ...);
    for (int j = 0; j < (int)_numNeighbours; ++j) {
        _outNeighbourIndices(i, j) = neighbourIndices[j];
        _outNeighbourSquaredDistances(i, j) = neighbourSquaredDistances[j];
    }
}
```

Two mechanical changes only:
1. Local buffers move inside the loop body — each OpenMP thread gets its own
   copy on its stack. No synchronization needed.
2. Loop variable `i` changes from `unsigned int` to `int` — OpenMP parallel for
   requires a signed integer loop variable (C++ OpenMP spec §2.9.2).

Output writes `_outNeighbourIndices(i, j)` and `_outNeighbourSquaredDistances(i, j)`
are indexed by row `i` — distinct rows per thread, no write conflict.

The KD-tree (`_kdTree->index->findNeighbors`) is read-only during the query
phase and nanoflann's KD-tree is documented as safe for concurrent read-only
queries.

The KD-tree build (`set_source_points` → `buildIndex`) remains serial. It is
called once per ICP iteration before the query loop and is not the measured
bottleneck.

## Build Integration

OpenMP is not currently linked in the project. Add to `library/CMakeLists.txt`:

```cmake
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(meshmonk_lib PUBLIC OpenMP::OpenMP_CXX)
endif()
```

`PUBLIC` scope propagates transitively to `_meshmonk_core`. The `if` guard
keeps the build working on macOS Apple Clang (which ships without OpenMP by
default) — the code path is unchanged when OpenMP is absent, just serial.

On macOS, libomp can be installed via `brew install libomp` and the build
configured with:
```
CMAKE_C_COMPILER=gcc-14 CMAKE_CXX_COMPILER=g++-14
```
or by passing `-DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp" -DOpenMP_CXX_LIB_NAMES=omp`.

## Profiling Driver

The profiling driver pins `OMP_NUM_THREADS=1` to match the ADR-002 baseline
convention. This pin means the parallelized build produces identical profiling
output to the serial build when run via the driver — **profiling a parallelized
build with OMP_NUM_THREADS=1 is safe and expected**. A separate profiling run
with `OMP_NUM_THREADS` unset (or set to the hardware thread count) is needed to
measure the actual parallel speedup.

## Validation

### Step 1 — Correctness: existing test suite (bitwise identical)

Because the parallelized loop has no floating-point reductions — each thread
writes to distinct rows of `_outNeighbourIndices` and
`_outNeighbourSquaredDistances` and reads only from the shared read-only
KD-tree — results must be **bitwise identical** to the serial build, not just
within tolerance.

Run the standard test suite (excluding golden, which is nightly-only):

```
pytest tests/ -q --ignore=tests/test_golden.py
```

A single failing test is a correctness regression. There is no tolerance budget
here: any divergence means the parallelization introduced a bug (most likely a
shared-buffer race from a buffer accidentally left outside the loop).

### Step 2 — Thread safety: ThreadSanitizer

Build with TSAN and run the fast integration tests:

```
cmake -DCMAKE_CXX_FLAGS="-fsanitize=thread -g -O1" ...
pytest tests/test_registration_e2e.py tests/test_e2e.py -q
```

TSAN reports data races at runtime. Any TSAN report is a blocker — a race that
doesn't manifest in correctness testing on a given run can manifest under load
or on a different scheduler. Common false positive: Eigen's internal thread-local
state; those can be suppressed with a TSAN suppression file if confirmed
benign, but should be investigated first.

TSAN and OpenMP interact: use `libomp` (LLVM) rather than `libgomp` (GCC) for
TSAN builds, or GCC with `-fsanitize=thread` and `OMP_NUM_THREADS=N` explicitly
set, since TSAN doesn't fully understand GCC's OpenMP runtime internals on all
platforms.

### Step 3 — Performance: profiling rerun

Run the profiling driver with `OMP_NUM_THREADS` set to the hardware thread count
(unpin from the driver's default of 1):

```
OMP_NUM_THREADS=$(nproc) python -m profiling.run_profile \
  --tiers 1k,10k,100k \
  --modes nonrigid,pyramid \
  --runs 5 --warmup 1 --seed 42 \
  --out docs/perf/hotspot-profile-YYYYMMDD-<sha>-parallel.md
```

Compare the new `NeighbourFinder::update` ms/call and share against
`hotspot-profile-20260423-0e3072f.md`. The correspondence labels should shrink
proportionally to the thread count (up to the Amdahl ceiling). If they do not
shrink — or share actually stays the same — the OpenMP linkage is likely missing
(the `#pragma` was silently ignored).

To confirm OpenMP is active at runtime, add a one-liner check to the driver or
verify via:
```python
import meshmonk
# or from C++ side: omp_get_max_threads() > 1
```

### Step 4 — Performance gate: benchmark harness

Run the benchmark harness against `baseline.json` under the ADR-002 perf-track
protocol:

```
OMP_NUM_THREADS=$(nproc) python benchmarks/bench_nonrigid.py
OMP_NUM_THREADS=$(nproc) python benchmarks/bench_pyramid.py
```

Gate: ≥5% wall-clock improvement at one or more mesh sizes compared to
`baseline.json`. The baseline was captured without OpenMP (serial build), so any
measurable parallel speedup should exceed 5% at 10K+ vertices given the
correspondence labels' combined 50–65% share.

If the gate is not met — possible if the test machine has only 2 physical cores
or if fork-join overhead dominates at 1K — document the result in the bead
notes and reassess the schedule strategy (D3 in ADR-005).

### Step 5 — Update baseline

If the benchmark gate passes, update `baseline.json` to capture the new
parallel-build wall-clock as the reference for future perf beads:

```
python benchmarks/bench_nonrigid.py --update-baseline
python benchmarks/bench_pyramid.py --update-baseline
```

The new baseline represents a parallelized build. All future perf beads must
build with OpenMP linked to produce comparable numbers.
