# ADR-007: NeighbourFinder Warm-Start Cache (Bead 729)

**Status:** Accepted
**Date:** 2026-04-24
**Parent:** [ADR-006 NeighbourFinder profiling](./ADR-006-neighbourfinder-profiling.md)
**Decision tree:** [bead 5jf follow-up](#decision-tree-context)
**Evidence:** before [hotspot-profile-20260424-(baseline)](../../docs/perf/neighbourfinder-breakdown-20260423-239b4cb.md); after [neighbourfinder-caching-20260424.md](../../docs/perf/neighbourfinder-caching-20260424.md)

## Outcome

`NeighbourFinder::update()` now pre-seeds nanoflann's `KNNResultSet::worstDist`
with an inflated ceiling derived from the previous iteration's neighbour
indices (distances recomputed against current positions). The tighter initial
`worstDist` prunes sibling branches earlier in the kd-tree traversal,
reducing average `tree_query` cost at 100K nonrigid from **~200 ns → ~100 ns
per query (≈50% reduction)**. Per-iteration `NeighbourFinder::update` drops
**~19%** and full nonrigid wall-time drops **~6%**, with bit-identical output
to the uncached path (verified by `test_cached_matches_uncached_nonrigid`).

## Context

ADR-006 established that `tree_query` (nanoflann `findNeighbors`) is the
dominant cost bucket inside `NeighbourFinder::update()`. The 5jf decision
tree recommended neighbour caching across ICP iterations as the next
investigation: at convergence, the query mesh moves only slightly between
iterations, so the previous iteration's k-NN remain close to optimal. An
upper bound on the true k-th distance — i.e., the *max* of cached-index
distances evaluated at current positions — is free correctness information
that lets the tree prune more aggressively.

## Decisions

### D1: Ceiling-only warm-start (not result-seeding)

**Firmness: FIRM**

The cache tightens `worstDist` by overwriting the sentinel slot
`neighbourSquaredDistances[k-1]` directly after `KNNResultSet::init()`. The
cache does **not** pre-populate the result set with cached (index, distance)
entries.

**Rationale:** `nanoflann::KNNResultSet::addPoint` has no duplicate-index
rejection. If cached entries are pre-added, the tree search re-visits those
same source indices, computes the same distances, and calls `addPoint` again
→ duplicate indices land in the result. Ceiling-only avoids this entirely:
no entries are added from the cache, and the tree search naturally fills
the result with the true k-NN (cached or not).

**Correctness argument (proof sketch):** Let `S` be the cached k indices,
`T` the true k-NN at current positions. Any point not in `T` has distance
`≥ t_k` (the true k-th). If `S ≠ T`, then `S \ T` is non-empty, and every
index in `S \ T` has current distance `≥ t_k`. Therefore
`c_k = max_i dist(query, source[S_i]) ≥ t_k`. The ceiling is always a valid
upper bound on the true k-th distance, so no true k-NN branch is pruned.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **Ceiling-only (chosen)** | No duplicates; minimal code; correct-by-construction | Cached point with `dist == ceiling` would fail nanoflann's strict-`<` check — handled by D2 |
| Pre-seed result entries | Intuitive | Creates duplicate indices when tree re-visits (nanoflann's `addPoint` has no uniqueness check) |
| Custom `KNNResultSet` with dedup | No duplicates | Header-only nanoflann is vendored as upstream source — patching it risks diverging from upstream |

### D2: Inflate ceiling by `(1 + eps)` to cancel nanoflann's approximate-search multiplier

**Firmness: FIRM**

The ceiling is written as `c_k * (1.0f + searchEps)` where `searchEps`
matches the `eps` passed to `SearchParams(32, searchEps, true)`.

**Rationale:** nanoflann's branch-descent decision at `KDTreeSingleIndexAdaptor::searchLevel`
(vendor/nanoflann.hpp:1399) is

```cpp
if (mindistsq * epsError <= result_set.worstDist())  descend;
```

with `epsError = 1 + eps = 1.0001`. With `worstDist = c_k` (raw), branches
are pruned when `mindistsq > c_k / 1.0001 ≈ c_k * 0.9999` — strictly tighter
than the exact boundary `mindistsq > c_k`. Because `mindistsq` is only a
*lower bound* on distances in a branch, a branch containing a true k-NN
point with dist `∈ (c_k * 0.9999, c_k]` could be pruned.

Inflating the ceiling to `c_k * (1 + eps)` cancels the multiplier: the prune
condition becomes `mindistsq * 1.0001 > c_k * 1.0001` ⟺ `mindistsq > c_k`,
the exact boundary. This restores bit-identical output across cached and
uncached paths.

**What this cost before fixing:** With the raw-ceiling variant, nonrigid
registration at 100K accumulated a max vertex-position drift of ~5 units
over 5 ICP iterations versus the uncached path — verified by comparing
`registration_result` against the stale regression reference (max drift
5.15 vs. 2.85e-4 for the uncached build).

### D3: Cache invalidation on source-pointer change

**Firmness: FIRM**

`NeighbourFinder::set_source_points()` invalidates the cache iff the source
pointer or row count changes. Same pointer + same row count → keep cache.

**Rationale:** Cached indices are row-indices into the source matrix. The
repo's call pattern (push filter: static target mesh; pull filter: floating
mesh mutated in place via persistent pointer) preserves pointer identity
across ICP iterations, so the cache survives. A different pointer could
mean a completely different mesh and we must be conservative — drop the
cache. A changed row count is a definitive invalidation.

**What this does not handle:** `set_parameters()` with a new `k` invalidates
the cache separately (row shape of `_cachedNeighbourIndices` is
`numQueried × k`). A call-site that reuses a `NeighbourFinder` instance
across two unrelated meshes with identical row count + pointer reuse would
keep stale cache — not a pattern used in this codebase, and defensive
validation is not warranted on hot paths.

### D4: Test-only Python toggle `_set_neighbour_caching`

**Firmness: FIRM**

A process-global default `registration::g_neighbour_caching_default` (true
by default) is exposed via the underscore-prefixed bindings
`meshmonk._set_neighbour_caching` and `meshmonk._get_neighbour_caching`.
Each `NeighbourFinder` copies this default at construction, so the toggle
must be set **before** invoking register functions.

**Rationale:** The correctness criterion ("cached and uncached paths
produce identical output") is only meaningfully testable with a runtime
toggle. Underscore-prefix signals this is test infrastructure, not a stable
public API. Keeping the flag module-global (not a kwarg on register) avoids
cluttering the public signatures for a debug hook.

## Decision tree context

From ADR-006 / bead 5jf:

> **Recommended next bead:** neighbour caching across ICP iterations
> (#1 in the 5jf decision tree). The `NeighbourFinder::tree_query`
> ScopedTimer retained in the code provides a direct before/after signal
> for that investigation.

This ADR is the realization of that recommendation. The retained
`tree_query` label did provide the signal — per-call cost dropped below
the 5× noise floor (now flagged UNRELIABLE), which is the strongest
observable form of "faster" at this granularity.

## Follow-ups

- **Batch SIMD queries (#2 in the 5jf decision tree)** remains deferred.
  Tree traversal is now fast enough that the next bucket to investigate is
  `SymmetricCorrespondenceFilter::update` post-NF work (triplet-list
  assembly + sparse matrix construction in `_update_affinity` / `fuse_affinities`),
  which now dominates after NF dropped.
- **Skip redundant tree rebuilds for push filter** (target is static
  across ICP): separate optimization not covered here. `NeighbourFinder::set_source_points`
  unconditionally rebuilds the kd-tree even when pointer + data are
  identical. Bead-sized opportunity if re-examined.
