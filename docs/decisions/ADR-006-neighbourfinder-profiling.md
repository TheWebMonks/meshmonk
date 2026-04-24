# ADR-006: NeighbourFinder Line-Level Profiling

**Status:** Accepted
**Date:** 2026-04-23
**Design:** [NeighbourFinder profiling design doc](../../history/2026-04-23-neighbourfinder-profiling-design.md)
**Parent:** [ADR-004 profiling infrastructure](./ADR-004-profiling.md)
**Related:** [ADR-005 NeighbourFinder OpenMP](./ADR-005-neighbourfinder-openmp.md) (superseded — not implemented)
**Evidence:** [NeighbourFinder breakdown 20260423](../../docs/perf/neighbourfinder-breakdown-20260423-239b4cb.md)

## Outcome

`tree_query` (`_kdTree->index->findNeighbors`) is the dominant cost inside
`NeighbourFinder::update()`. At 100K nonrigid it accounts for ~15% of
registration wall-time (the `NeighbourFinder::update` share); all other
per-vertex sub-buckets (`buffer_alloc`, `query_setup`, `result_copy`) are
below the 5× noise floor (UNRELIABLE per ADR-004 D2) and have been removed
from instrumentation.

**Recommended next bead:** neighbour caching across ICP iterations (#1 in
the 5jf decision tree). The `NeighbourFinder::tree_query` ScopedTimer
retained in the code provides a direct before/after signal for that
investigation.

---

## Context

ADR-005 (OpenMP parallelization) was closed without implementation — prior
bead `9f5` measured a regression (0.93–0.99×) at 7K. Before committing to
any further optimization, we needed to know WHERE inside
`NeighbourFinder::update()` the cost lives. The four candidate buckets were:
`buffer_alloc`, `query_setup`, `tree_query`, `result_copy`.

This bead added per-bucket `ScopedTimer` instrumentation, ran at 100K
nonrigid (5 measured repetitions, 1 warmup, median across reps), and
established that `tree_query` dominates while the other three buckets are
below the noise floor.

---

## Decisions

### D1: RAII scope bracketing for each cost bucket

**Firmness: FIRM**

Each instrumented cost bucket uses a nested `{}` scope to trigger RAII
destruction at the right point. Explicit destructor calls (`_t.~ScopedTimer()`)
are not used — they produce undefined behavior after the destructor runs.

**Rationale:** Correct C++. The scoped pattern matches existing usage
throughout the codebase (ADR-004 D2 naming convention). Investigation
complete: three of the four scopes have now been removed (see D4); only
`tree_query` remains.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **RAII scope brackets (chosen)** | Correct C++, matches existing pattern | Requires hoisting declarations above scope boundaries |
| Explicit destructor call | No variable hoisting needed | UB after first call; reader confusion |
| Timestamp pairs | No object lifetime concern | More boilerplate; doesn't match existing infrastructure |

**What would invalidate this:** A change to the `ScopedTimer` API that
provides a manual `.stop()` method with correct semantics. Until then, scope
brackets are the right tool.

---

### D2: Hoist buffer declarations above the loop

**Firmness: FIRM**

`queriedFeature`, `neighbourIndices`, `neighbourSquaredDistances`, and
`knnResultSet` are declared before the `for` loop. `resize()` and
reinitialization happen inside the loop on each iteration. After iteration 0,
`resize()` is a no-op (vectors stay at target capacity).

**Rationale:** Separates first-allocation cost (amortized over the loop) from
per-iteration work (what the scope timers measure). Without hoisting, the
`buffer_alloc` scope would capture both heap allocation (first iteration) and
a capacity check (all subsequent iterations) — making the first-rep
measurement an outlier and per-iteration numbers misleading. The hoist is a
positive structural change regardless of instrumentation: it documents that
the allocations are one-time.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **Hoist above loop (chosen)** | Clean first-alloc vs. per-iter separation | Declarations further from use |
| In-loop declaration | Declarations near use | Conflates first-alloc cost with per-iter cost in profiling |
| `thread_local` pre-sized buffers | Zero alloc every call | `KNNResultSet` must be reinit'd each query anyway; added lifecycle complexity |

**What would invalidate this:** A future refactor that moves the loop body
into a lambda or separate function, at which point the hoist may move back
inside the function scope naturally.

---

### D3: Run at 100K nonrigid only

**Firmness: FIRM**

The profiling investigation runs one nonrigid registration at 100K with
≥5 measured repetitions (1 warmup discarded), median across reps. No
multi-tier scaling fit is attempted for this investigation.

**Rationale:** The goal is to identify which bucket dominates — visible at the
tier where `NeighbourFinder` is measurably hot. Multi-tier runs are needed for
scaling-exponent fits, but the sub-bucket investigation has a single binary
question (which bucket?), not a scaling question. Running all three tiers
would triple the investigation time with no new information. The 100K tier
produces >5 million per-bucket invocations, well above the noise floor.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **100K nonrigid only (chosen)** | Fast; sufficient invocations for noise floor | No scaling exponent for sub-buckets |
| All three tiers (1K/10K/100K) | Full scaling picture | 3× longer; scaling question is not the purpose of this investigation |
| 10K only | Faster | Fewer invocations; sub-bucket shares may differ from production scale |

**What would invalidate this:** A follow-up investigation that specifically
asks HOW sub-bucket costs scale with n (e.g., to distinguish O(n) vs O(n log n)
allocation patterns). At that point, re-enable multi-tier runs.

---

### D4: Remove per-vertex sub-function ScopedTimer labels after investigation

**Firmness: FIRM**

Sub-function `ScopedTimer` labels inside per-vertex loops (e.g.
`NeighbourFinder::buffer_alloc`, `NeighbourFinder::query_setup`,
`NeighbourFinder::result_copy`) must be removed once the investigation that
added them is complete, unless they are gated behind a separate verbose build
flag (e.g. `-DMESHMONK_PROFILING_VERBOSE=ON`).

**Rationale:** The `MESHMONK_PROFILING=ON` build is a shared artifact used
by all profiling investigations — not a per-bead sandbox. Per-vertex
instrumentation inside loops with 100K+ iterations adds significant timer
overhead:

- 4 sub-labels × 100K vertices × 403 `NeighbourFinder` calls ≈ 161.2M
  `ScopedTimer` constructions per run.
- At 32–35 ns per construction, that is ~5.1–5.6 s of pure overhead per run.
- This contaminates the `NeighbourFinder::update` share measurement and
  inflates every label that nests inside it (e.g.
  `SymmetricCorrespondenceFilter::update`).

After removing the three noise-floor scopes and retaining only `tree_query`,
the overhead drops to ~40.3M constructions (~1.4 s), isolated to the one
label that is actually measurable.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **Remove after investigation (chosen)** | Zero ongoing overhead; clean shared build | Instrumentation must be re-added if the question recurs |
| Keep all four sub-labels permanently | No re-add cost | ~5s overhead per run contaminates every future investigation |
| Gate behind `-DMESHMONK_PROFILING_VERBOSE=ON` | Keeps code; opt-in overhead | Added CMake complexity; not yet needed — remove first, add gate only if question recurs |

**What would invalidate this:** A future investigation that again needs
sub-bucket breakdown inside a per-vertex loop. At that point, gate the labels
behind a `-DMESHMONK_PROFILING_VERBOSE=ON` flag rather than adding them
unconditionally.
