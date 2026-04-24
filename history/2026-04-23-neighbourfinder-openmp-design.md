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

## Review Pass 1 Findings

### ADRs reviewed

- **ADR-001 meshmonk-modernization** — strategic parent; no direct conflicts with this design but relevant for (a) D5 CMake + scikit-build-core build backend (constrains how OpenMP is wired in for wheels), (b) D2 nanobind/ABI strategy (parallel macOS wheels must bundle libomp — not currently handled by cibuildwheel config).
- **ADR-002 perf-track** — VERY RELEVANT. D1 benchmark-first, D2 OpenMP toolchain, D3 tolerance vs bitwise, D4 merge gate, D6 baseline refresh policy. Contains a **D4 Outcome section for bead 9f5 (2026-04-22)** that documents this exact optimization (OpenMP on `NeighbourFinder::update`) being **attempted, measured, and closed below gate**. The design and ADR-005 do not reference this prior attempt.
- **ADR-003 api-polish** — unrelated surface polish; no bearing here.
- **ADR-004 profiling** — provides the measurement infrastructure this design relies on; `run_profile.py` hard-pins `OMP_NUM_THREADS=1` at module-load time (ADR-004 D1 does not require this, but the driver enforces it — see Finding 3 below).
- **ADR-005 neighbourfinder-openmp** — paired ADR. Its D1–D5 are FIRM/FLEXIBLE as stated; specific conflicts listed in findings below.

### Findings

#### F1 [ARCH] CRITICAL: Prior-art gap — bead 9f5 already tried this exact change and closed below gate

ADR-002 lines 168–205 (the `D4 Outcome (2026-04-22): 9f5 closed below recalibrated gate` block) documents that `#pragma omp parallel for` was added to this exact loop in `NeighbourFinder::update()` on 2026-04-22, measured against the 1K/3K/7K benchmark harness with `OMP_NUM_THREADS` set to 4 and 14, and **regressed every scenario** (0.93–0.99× across nonrigid/pyramid/rigid @ 1K/3K/7K). The bead was closed, the code stashed, and the closure rationale was: fork-join overhead exceeds compute savings at the benchmark-reachable mesh sizes (per-thread workload ≤1750 pts at 7K ÷ 4 threads).

Neither the paired ADR-005 nor this design doc reference 9f5. The profiling-driven evidence in `hotspot-profile-20260423-0e3072f.md` differs from the 9f5 picture (it measures 1K/10K/100K via the profiling driver rather than 1K/3K/7K via the benchmark harness, and it shows super-linear scaling k≈1.26), which is a *new* reason to look at 100K — but 100K is not reachable from the benchmark harness (which caps at 7K, per `benchmarks/README.md`).

**Required change:** The design must explicitly cite 9f5's closure, explain what's different this time (only 100K-via-profiling-driver adds new information; 7K-via-harness is unchanged), and acknowledge that the "≥5% at 10K+ vertices" gate in Step 4 cannot be measured by the current benchmark harness. Also explicit: ADR-005 should add a "what would invalidate this" clause and/or an explicit "supersedes 9f5" reference so the history is discoverable.

#### F2 [ARCH] Amdahl ceiling misstated — the three correspondence labels are nested, not additive

The "Why Here" table claims the correspondence pipeline is ~50–65% of nonrigid wall-time and implies parallelizing NeighbourFinder lifts all three. This is a misreading of the profile — the labels are **nested**: `CorrespondenceFilter::update` *contains* `NeighbourFinder::update`; `SymmetricCorrespondenceFilter::update` *contains* two `CorrespondenceFilter::update` calls. The hotspot report states this explicitly (NeighbourFinder share=14.7% at nonrigid@100K, Amdahl ceiling **1.2×**), and the report's own candidate list notes the 1.2× ceiling per label.

The design's 50–65% number is the sum of the three labels, which double- and triple-counts NeighbourFinder's time. The attainable speedup for the nonrigid path from parallelizing NeighbourFinder is bounded by the **NeighbourFinder share alone** (≤14.7% at 100K → ceiling 1.17×), plus the small affinity-construction overhead (~10% of CorrespondenceFilter, ~1.6% of nonrigid wall-time) which is NOT parallelized by this change.

**Required change:** Rewrite "Why Here" to cite the raw NeighbourFinder share (6.7% pyramid@1K / 12.6% nonrigid@10K / 14.7% nonrigid@100K, Amdahl ceiling 1.1× / 1.1× / 1.2×) and set realistic expectations: on nonrigid@100K the best case is ~17% wall-clock improvement, not proportional-to-thread-count. The 9f5 closure already demonstrated that 14-thread fork-join cost per invocation (15–30 µs × ~400 invocations = 6–12 ms of scheduler overhead) consumes a meaningful fraction of the attainable ceiling at 7K; at 100K the invocations are still ~400 but per-invocation work is ~56 ms, so overhead is negligible there.

#### F3 CRITICAL: Step 3 profiling rerun command cannot actually unpin OMP_NUM_THREADS

`profiling/run_profile.py` lines 21–23 set `os.environ["OMP_NUM_THREADS"] = "1"` and `os.environ["EIGEN_DONT_PARALLELIZE"] = "1"` unconditionally at module load, before any other imports. The command `OMP_NUM_THREADS=$(nproc) python -m profiling.run_profile …` in Step 3 of the Validation section is silently overridden — the driver will still run single-threaded. Every measured share in the rerun will match the baseline, which this design treats as evidence that "the `#pragma` was silently ignored" — a false negative.

**Required change:** Either (a) patch `run_profile.py` to respect a `MESHMONK_UNPIN_OMP=1` environment override (or a `--threads N` CLI flag), (b) document a separate parallel-profiling driver script, or (c) document an in-process override like calling `omp_set_num_threads()` via a small helper before registration runs. Option (a) is cleanest given ADR-004 D1's "driver controls profiling lifecycle" stance.

#### F4 CRITICAL: Step 4 benchmark commands are wrong — benches are pytest files, baseline captured via conftest fixture

- `python benchmarks/bench_nonrigid.py` is not a runnable script (no `__main__`, no argparse) — it's a pytest module.
- `python benchmarks/bench_pyramid.py --update-baseline` — no such flag exists anywhere.
- `benchmarks/conftest.py` has a **session-scoped autouse `pin_threads` fixture** that sets `OMP_NUM_THREADS=1` for the entire run, so `OMP_NUM_THREADS=$(nproc) uv run pytest benchmarks/ …` will have `OMP_NUM_THREADS` overwritten to 1 at session start. A parallel benchmark run needs either the fixture disabled/gated behind an env var (e.g. `MESHMONK_BENCH_UNPIN=1`) or a separate `benchmarks-parallel/` directory.
- Baseline refresh uses `--benchmark-json=benchmarks/baseline.json`, not `--update-baseline` (see `benchmarks/README.md` lines 17–19).
- No benchmark tier exists above 7K; `benchmarks/README.md` cap-at-7K decision (lines 33–67) is explicit. The design's "≥5% at 10K+ vertices" criterion is unreachable at the benchmark harness.

**Required change:** Replace Steps 4–5 with working invocations against the actual pytest harness AND either (a) add an unpinning mechanism to conftest, or (b) acknowledge that the parallel-path benchmark must use the profiling driver's 10K/100K meshes instead, not the benchmark harness. Given ADR-002 D4 is the documented merge gate, clarify which harness satisfies the gate for this bead. The "≥5% at 10K+ vertices" threshold in Step 4 must be relocated to the profiling driver comparison (Step 3) because the benchmark harness does not reach 10K.

#### F5 [ARCH] CONFLICTS WITH ADR-002 D2: `find_package(OpenMP)` without REQUIRED

ADR-002 D2 (FIRM) prescribes `find_package(OpenMP REQUIRED)`. This design's Build Integration section and ADR-005 D4 (FIRM) prescribe `find_package(OpenMP)` without REQUIRED, so the build "silently" falls back to serial on macOS without libomp. This is a deliberate softening of an ADR-002 FIRM decision, not a silent deviation — but it must be flagged.

ADR-002 D2's "What would invalidate this" clause does note that macOS libomp availability was an open question at write time ("verify during 9f5"). The pragmatic softening (guarded find_package) is reasonable but needs an explicit "this amends ADR-002 D2" note in ADR-005, or an update to ADR-002 D2. Without either, the two ADRs read inconsistently.

**Required change:** ADR-005 D4 must add: "This amends ADR-002 D2 (originally `REQUIRED`) to `find_package(OpenMP)` guarded, to preserve macOS Apple Clang default-build." ADR-002 D2 should get a cross-reference noting the amendment.

#### F6 [ARCH] CONFLICTS WITH ADR-002 D3: "bitwise identical" claim is wrong and conflicts with the stated invariance policy

Step 1 of Validation claims that the parallelized loop produces "bitwise identical" output because there are no floating-point reductions. This is narrowly true for THIS loop (each thread writes distinct rows; no reduction). However:

- ADR-002 D3 (FIRM) explicitly prescribes **tolerance-based, not bitwise** invariance for OpenMP changes (D3 rationale: "9f5 and bdt … may produce output that differs by a few ULPs on some entries"). ADR-002 D4 criterion 5 scopes bitwise only to 6a5 (memory layout, no new operations).
- The existing e2e test suite (`tests/test_registration_e2e.py`) uses scalar aggregate thresholds (`step.mean() < 10.0` etc., per ADR-002 D3 paragraph 2), NOT per-vertex equality. A bitwise regression would not be caught by running this suite; the suite "passing" is neither necessary nor sufficient evidence of bitwise identity.
- The actual bitwise guard in the repo is `tests/golden/test_legacy_baseline.py` and `tests/test_memory_layout_regression.py` (per ADR-002 D4 Outcome 6a5). Neither is mentioned in Step 1.
- The 9f5 closure explicitly demonstrated bitwise-identical aligned vertex coordinates under OMP=1 vs OMP=4 (max-abs-diff = 0.0), so the claim that *this specific loop* preserves bits is empirically supported — but the mechanism for *verifying* it in Step 1 is wrong.

**Required change:** Either (a) drop the "bitwise identical" framing and align with ADR-002 D3 tolerance-based language, with the golden-path test as the tight bound; or (b) keep the bitwise claim but point explicitly to `tests/test_memory_layout_regression.py` + `tests/golden/test_legacy_baseline.py` as the verification mechanism, and add a new test specifically comparing OMP=1 vs OMP=N aligned-vertex arrays (parity test per ADR-002 D3 paragraph 2), which is the guard 9f5 actually used.

#### F7 [ARCH] Column-major output matrix + row-striped parallelism = false sharing risk

`MatDynInt` and `MatDynFloat` are declared `Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>` in `library/src/NeighbourFinder.hpp` lines 9–11. Eigen's default storage order for `Dynamic×Dynamic` matrices is **column-major**. With `schedule(static)`, adjacent iterations `i`, `i+1` (same column `j`) map to adjacent memory addresses at stride 1 float (or 1 int). Two threads processing adjacent chunks write to cache lines that straddle their boundary → classic false sharing on the last ~16 floats of thread T's chunk and the first ~16 of thread T+1's chunk.

At k=3 neighbours and `_numQueriedElements = 100K`, each column is ~400 KB contiguous. At static scheduling with 14 threads → ~7K rows/thread → each thread owns a ~28 KB run in each column. False-sharing pressure is confined to the chunk boundaries (7K rows in × 4 bytes = 28 KB chunks; cache line is 64 B; only 14 boundaries × 3 columns = 42 contested lines per pass). This is small relative to the compute, so likely not measurable — but the design claim "no write conflict" (line 96) is incomplete without noting the column-major layout and the false-sharing characterization.

**Required change:** Add a sentence in "The Change" acknowledging that output matrices are column-major; static-scheduled row partitioning yields contiguous per-thread regions so false sharing is confined to chunk boundaries and is expected to be negligible. Alternatively, if future measurement shows it, the output matrices could be declared `Eigen::RowMajor` — at the cost of an Eigen API audit for all consumers (`get_indices()`, `get_distances()`, Python bindings).

#### F8 Hotspot share table in "Why Here" cites percentages not present in the cited profile

The design's "Why Here" table shows:
- `SymmetricCorrespondenceFilter::update` 22%
- `CorrespondenceFilter::update` 16%
- `NeighbourFinder::update` 15%

The hotspot profile `hotspot-profile-20260423-0e3072f.md` at nonrigid@100K shows 22.1%, 16.0%, 14.7% respectively. Small rounding, acceptable. But the table header says "Share (nonrigid, 100K)" so the 100K picture is in scope — and at 100K the ADR-002 D4 Outcome for 9f5 has a deferred-revisit clause: "at 50K points split across 14 threads ≈ 3600 pts/thread, compute-per-invocation should clear the fork-join break-even. The gate to clear at that point is the original `≥2×` target". This design's Step 4 gate is ≥5%, inconsistent with ADR-002 D4's deferred-revisit stipulation of ≥2×.

**Required change:** Reconcile the Step 4 merge gate with ADR-002 D4's stipulated re-entry threshold (original `≥2×` on 5K-nonrigid, deferred to 50K-exists). The 5% threshold in this design is less stringent than ADR-002 specifies; either this design must explicitly amend ADR-002 D4 (with rationale — e.g. the profiling driver's 100K hotspot evidence justifies a softer gate given Amdahl 1.2×), or adopt the ≥2× target. Note that given F2's correct Amdahl ceiling of ~1.2× for the NeighbourFinder share alone, a ≥2× gate is arithmetically unreachable and the design should argue for a new gate with the Amdahl bound as the ceiling.

#### F9 TSAN + libgomp guidance is understated

Step 2 says "TSAN doesn't fully understand GCC's OpenMP runtime internals on all platforms" — understated. The official TSAN + libgomp interaction is: libgomp does not provide TSAN annotations, so **every OpenMP barrier and task boundary reports false-positive races**. The practical workarounds are:
1. Build libgomp with TSAN annotations (`--disable-linux-futex` + source patches, non-trivial).
2. Use libomp (LLVM) which has TSAN annotations.
3. Use `archer` (LLVM's OpenMP-aware TSAN wrapper).

The "GCC with `-fsanitize=thread` and `OMP_NUM_THREADS=N` explicitly set" workaround does not fix the false positives; it only makes thread count deterministic. The libomp recommendation is the only one that actually works.

**Required change:** Remove the GCC-with-`-fsanitize=thread` workaround as an alternative; keep only the libomp path. Link to `archer` as an optional enhancement.

#### F10 Runtime OpenMP detection path is loose

Step 3's "To confirm OpenMP is active at runtime" section is incomplete:
```python
import meshmonk
# or from C++ side: omp_get_max_threads() > 1
```
There is no Python binding for `omp_get_max_threads()` (none of the five profiling bindings in ADR-004 D4 exposes it). A passing `OMP_NUM_THREADS=N python -c "..."` does not confirm the `_meshmonk_core.so` was linked against OpenMP — it only confirms the environment variable is set.

**Required change:** Add a concrete nanobind-exposed probe: e.g. `meshmonk.openmp_active() -> bool` (calls `#ifdef _OPENMP return true; #else return false; #endif`), returned alongside `omp_get_max_threads()` if enabled. This is a ~5-line binding addition and the only robust way to verify the compiled artifact has OpenMP linked. Without it, Step 3's "If [shares] do not shrink … the OpenMP linkage is likely missing" is speculation.

#### F11 Missing: Eigen OpenMP oversubscription risk

`meshmonk_lib` uses Eigen (vendored). Eigen itself uses OpenMP when linked — this is flagged in ADR-002 D2 ("Eigen already links OpenMP when available"). With OpenMP enabled project-wide and no `Eigen::setNbThreads(1)` binding exposed (per `benchmarks/conftest.py` lines 60–64, confirmed), any Eigen matrix operation inside the parallel region (or in surrounding code) can spawn nested OpenMP threads → oversubscription. The public `OMP_NUM_THREADS` thread knob (ADR-002 D5) doesn't distinguish meshmonk-level from Eigen-level.

While the body of `NeighbourFinder::update` doesn't call Eigen kernels (only element-access `_outNeighbourIndices(i,j)` which is trivial), the broader implication for the correspondence stack (which does touch Eigen via affinity construction) is unaddressed. ADR-002 D2 says this is acceptable ("Eigen already probes OpenMP"), so it's not necessarily a problem — but the design should explicitly state that `OMP_NESTED` is OFF by default, Eigen operations at outer levels see `OMP_NUM_THREADS` threads, and oversubscription is therefore controlled.

**Required change:** Add a "Oversubscription" subsection in Build Integration confirming: Eigen's opportunistic OpenMP uses the same `OMP_NUM_THREADS` budget, `OMP_NESTED=0` (default) prevents nested spawn, and the batch-workload caveat from `benchmarks/README.md` lines 108–118 carries forward.

#### F12 Missing edge case: `_numQueriedElements == 0`

`set_queried_points` followed by `update()` with a zero-row `_inQueriedPoints` produces an empty output (no rows to write). The parallelized version still enters the parallel region but the loop body doesn't execute (upper bound 0). Behaviorally identical, but the fork-join cost is paid — ~15 µs per the profile — where the serial version paid only the loop-overhead check. This is a micro-regression on degenerate inputs, not a correctness issue, but should be noted. Guarding with `if (_numQueriedElements == 0) return;` is a one-liner.

**Required change:** Add to "The Change": handle `_numQueriedElements == 0` explicitly (early return) so the degenerate case doesn't pay fork-join cost.

#### F13 Minor: `int` loop variable truncation

The loop variable change from `unsigned int i` to `int i` is correct per OpenMP C++ requirements, but `(int)_numQueriedElements` truncates for meshes ≥ 2^31 vertices. Unreachable in practice (2.1B vertices ≈ 50GB of positions alone) but worth a `static_assert` or a comment.

**Required change:** Add a `// static_assert on reasonable upper bound` comment, or use `long` / `ssize_t` instead of `int` (OpenMP 3.0+ permits signed integer types).

#### F14 Minor: `#pragma omp parallel for` and MSVC

MSVC's OpenMP is 2.0-level (no C++ loop collapse etc.); `parallel for` with `schedule(static)` works fine. MSVC does not treat unknown pragmas as errors by default — `/W4` warns, but the build still succeeds without OpenMP. No action needed; the design's guarded-link approach already handles this. Just calling out that the paragraph "compilers treat unknown pragmas as warnings, not errors" is accurate for MSVC, Clang, GCC default settings.

#### F15 Minor: macOS wheel distribution (cibuildwheel + libomp)

ADR-002 D2 says "cibuildwheel already handles [macOS libomp]", but there's no visible cibuildwheel config in the repo that actually installs libomp for the macOS wheel build step. If the CI wheel build runs on macOS and doesn't have libomp, the guarded `find_package(OpenMP)` will silently build serial wheels. Users installing on Linux get parallel; users installing on macOS get serial. This is ADR-005 D4's explicit rationale, but the wheel-distribution side-effect isn't highlighted.

**Required change:** Add a sentence to Build Integration / ADR-005 D4 "What would invalidate this" clause: macOS wheels ship serial by default; if that's unacceptable at some future release, add libomp installation to the macOS cibuildwheel config and bump the macOS build step.

#### F16 Minor: `schedule(static)` assumption of uniform per-vertex cost

ADR-005 D3 is FLEXIBLE and acknowledges this. The assumption "k-NN query cost per vertex is approximately constant" is not strictly true — KD-tree `findNeighbors` prunes via `worstDist()`, and vertices in denser regions prune earlier (smaller search). In practice the variation is small for typical face meshes. Not a required change, just a note that ADR-005 D3's flexibility clause is the right escape hatch if measurement shows imbalance on a given mesh.

### Verdict

**BLOCKED: [F1, F3, F4 must be resolved before implementation]**

The design as written will not produce a valid measurement. F1 (9f5 prior art) is the load-bearing issue: this is the same change that was closed below gate 1 day before this design was written, and the design doesn't cite the prior attempt or explain what's different. F3 (profiling driver pins OMP=1) and F4 (benchmark harness pins OMP=1, bench files not scripts, no 10K+ tier) mean Step 3 and Step 4 of the Validation plan cannot execute as written — every measurement will be single-threaded, producing false-negative results that the design would interpret as "OpenMP not linked".

Non-blocking but important: F2 (Amdahl math off by 3×), F5/F6/F8 (three ADR conflicts with ADR-002 that need explicit amendment), F7 (false-sharing characterization), F10 (need a real runtime detection probe). F9, F11–F16 are refinements.

Recommended resolution path:
1. Add a prior-art section citing 9f5 closure; explain why the profiling-driver 100K evidence justifies a second attempt; specifically limit the measurement scope to the profiling driver at 100K (not the 7K benchmark harness).
2. Fix the measurement plan: patch `run_profile.py` to respect an unpin-override, OR add a standalone parallel-profiling script, OR document that the bead includes an infrastructure change to both `conftest.py` and `run_profile.py`.
3. Re-align the Amdahl ceiling expectation (F2) and the merge gate (F8): ≥5% at 100K is plausible; ≥2× is arithmetically unreachable given the 14.7% share.
4. Amend ADR-005 to explicitly supersede 9f5's closure, soften ADR-002 D2 (FIRM → amend), and align with ADR-002 D3 on tolerance-based invariance (F6).
