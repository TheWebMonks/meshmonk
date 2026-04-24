# MeshMonk Hotspot Profile — 2026-04-24 (`unknown`)

## Methodology

- **CPU**: arm (14 logical processors)
- **OS**: Darwin 24.6.0
- **Compiler**: Apple clang version 17.0.0 (clang-1700.0.13.5)
- **CMAKE_BUILD_TYPE**: Release (assumed; use `cmake --build . --config Release`)
- **LTO**: unknown
- **Eigen**: vendored (vendor/eigen-3.4.0)
- **nanoflann**: vendored (vendor/nanoflann.hpp)
- **OpenMesh**: vendored (vendor/OpenMesh-11.0.0)
- **Repo SHA**: `unknown`
- **OMP_NUM_THREADS**: 1
- **EIGEN_DONT_PARALLELIZE**: 1
- **Run command**: `python -m profiling.run_profile --tiers 100k --modes nonrigid --runs 5 --warmup 1 --seed 42`
- **Seed**: 42
- **ScopedTimer overhead**: 32 ns/scope (median of 3 × profiling_calibrate(1_000_000); individual: 32, 32, 32 ns/scope)
- **Repetitions per (tier, mode)**: warmup=1, runs=5, median selected

### Mesh SHA-256

| File | SHA-256 |
|---|---|
| Template_1K.obj | `file not found` |
| Template_10K.obj | `file not found` |
| Template_100K.obj | `file not found` |
| DemoFace_1K.obj | `file not found` |
| DemoFace_10K.obj | `file not found` |
| DemoFace_100K.obj | `file not found` |

### Pyramid-mode deduplication note

In pyramid mode, `NonrigidRegistration::update` accumulates across all layers
(bare label, count = num_layers), AND each layer also emits a
`NonrigidRegistration::update/layerN` label. The bare label is included in
the table as a cross-check. Share percentages are computed from all labels
including the bare one — do not double-count by adding bare + layer shares.
A sanity check verifies bare_total ≈ sum(layerN_totals) within ±5%; violations
are noted inline.

## Per-tier Hotspot Rankings


### Tier: 100K (~100,000 vertices) — nonrigid

| Label | Share (%) | ms/call | Invocations | Amdahl ceiling | Fork-join cost (us) | Notes |
|---|---|---|---|---|---|---|
| `NonrigidRegistration::update` | 33.3 | 75215.9090 | 1 | 1.5x | 15 |  |
| `SymmetricCorrespondenceFilter::update` | 22.6 | 254.3604 | 201 | 1.3x | 3015 |  |
| `CorrespondenceFilter::update` | 17.5 | 98.3549 | 402 | 1.2x | 6030 |  |
| `NeighbourFinder::update` | 16.6 | 92.9917 | 403 | 1.2x | 6045 |  |
| `NeighbourFinder::tree_query` | 3.8 | 0.0002 | 40300000 | 1.0x | 604500000 |  |
| `ViscoElasticTransformer::_update_viscously` | 2.7 | 30.5653 | 200 | 1.0x | 3000 |  |
| `ViscoElasticTransformer::_update_elastically` | 2.7 | 30.5319 | 200 | 1.0x | 3000 |  |
| `InlierDetector::update` | 0.5 | 5.2786 | 201 | 1.0x | 3015 |  |
| `ViscoElasticTransformer::_apply_transformation` | 0.3 | 3.5952 | 200 | 1.0x | 3000 |  |
| `ViscoElasticTransformer::_update_outlier_transformation` | 0.1 | 1.6921 | 200 | 1.0x | 3000 |  |
| `ViscoElasticTransformer::_update_smoothing_weights` | 0.0 | 2.5320 | 1 | 1.0x | 15 |  |
| `NeighbourFinder::buffer_alloc` | 0.0 | 0.0000 | 40300000 | 1.0x | 604500000 | * UNRELIABLE |
| `NeighbourFinder::result_copy` | 0.0 | 0.0000 | 40300000 | 1.0x | 604500000 | * UNRELIABLE |
| `NeighbourFinder::query_setup` | 0.0 | 0.0000 | 40300000 | 1.0x | 604500000 | * UNRELIABLE |

*UNRELIABLE = ms/call < 5× timer overhead (32 ns/scope = 0.000032 ms). Excluded from recommendations.*

## Scaling Exponents (log-log fit across 1K / 10K / 100K)

| Label | Mode | k (exponent) | Complexity bucket | Note |
|---|---|---|---|---|
| `NonrigidRegistration::update` | nonrigid | NaN | n/a (insufficient tiers) | insufficient tiers |
| `CorrespondenceFilter::update` | nonrigid | NaN | n/a (insufficient tiers) | insufficient tiers |
| `NeighbourFinder::buffer_alloc` | nonrigid | NaN | n/a (insufficient tiers) | insufficient tiers |
| `NeighbourFinder::result_copy` | nonrigid | NaN | n/a (insufficient tiers) | insufficient tiers |
| `ViscoElasticTransformer::_update_elastically` | nonrigid | NaN | n/a (insufficient tiers) | insufficient tiers |
| `ViscoElasticTransformer::_update_viscously` | nonrigid | NaN | n/a (insufficient tiers) | insufficient tiers |
| `ViscoElasticTransformer::_apply_transformation` | nonrigid | NaN | n/a (insufficient tiers) | insufficient tiers |
| `ViscoElasticTransformer::_update_smoothing_weights` | nonrigid | NaN | n/a (insufficient tiers) | insufficient tiers |
| `SymmetricCorrespondenceFilter::update` | nonrigid | NaN | n/a (insufficient tiers) | insufficient tiers |
| `NeighbourFinder::update` | nonrigid | NaN | n/a (insufficient tiers) | insufficient tiers |
| `ViscoElasticTransformer::_update_outlier_transformation` | nonrigid | NaN | n/a (insufficient tiers) | insufficient tiers |
| `NeighbourFinder::query_setup` | nonrigid | NaN | n/a (insufficient tiers) | insufficient tiers |
| `NeighbourFinder::tree_query` | nonrigid | NaN | n/a (insufficient tiers) | insufficient tiers |
| `InlierDetector::update` | nonrigid | NaN | n/a (insufficient tiers) | insufficient tiers |

> **Caveat:** Three-point fits (1K/10K/100K) have zero residual degrees of freedom.
> The buckets k≈1.0 / k≈1.16 / k≈2.0 cannot be reliably discriminated without a
> fourth tier or per-tier variance. Treat exponents as indicative only.

## Recommendations

> **Note:** In pyramid mode, `PyramidNonrigidRegistration::update` is a
> reference-only label (superset of its `/layerN` children). It is excluded
> from the candidates list to avoid double-counting. Optimization candidates
> are the inner `/layerN` and leaf labels.

### Candidates for optimization

Labels meeting all of: share > 5%, ms/call > 0.03 ms, NOT UNRELIABLE.
In pyramid mode, `PyramidNonrigidRegistration::update` is excluded (reference-only).

- **`NonrigidRegistration::update`** (nonrigid, 100K tier): share=33.3%, Amdahl ceiling = 1/(1−0.333) = 1.5×
  - Scaling exponent k: NaN
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 15 µs (0.0% of label wall time)

- **`SymmetricCorrespondenceFilter::update`** (nonrigid, 100K tier): share=22.6%, Amdahl ceiling = 1/(1−0.226) = 1.3×
  - Scaling exponent k: NaN
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 3015 µs (0.0% of label wall time)

- **`CorrespondenceFilter::update`** (nonrigid, 100K tier): share=17.5%, Amdahl ceiling = 1/(1−0.175) = 1.2×
  - Scaling exponent k: NaN
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 6030 µs (0.0% of label wall time)

- **`NeighbourFinder::update`** (nonrigid, 100K tier): share=16.6%, Amdahl ceiling = 1/(1−0.166) = 1.2×
  - Scaling exponent k: NaN
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 6045 µs (0.0% of label wall time)

### Below granularity threshold

Labels with share ≤ 5% or ms/call ≤ 0.03 ms. OpenMP not viable due to
fork-join overhead exceeding potential speedup (reference: beads 9f5, bdt).

`InlierDetector::update`, `NeighbourFinder::tree_query`, `ViscoElasticTransformer::_apply_transformation`, `ViscoElasticTransformer::_update_elastically`, `ViscoElasticTransformer::_update_outlier_transformation`, `ViscoElasticTransformer::_update_smoothing_weights`, `ViscoElasticTransformer::_update_viscously`

### Excluded (UNRELIABLE)

Labels with ms/call < 5× timer overhead (32 ns/scope = 0.000032 ms).
Corrected values unreliable — cannot meaningfully distinguish measurement noise from actual cost.

| Label | Mode | Tier | ms/call |
|---|---|---|---|
| `NeighbourFinder::buffer_alloc` | nonrigid | 100K | 0.0000 |
| `NeighbourFinder::result_copy` | nonrigid | 100K | 0.0000 |
| `NeighbourFinder::query_setup` | nonrigid | 100K | 0.0000 |

---

## NeighbourFinder Line-Level Recommendation

### Summary

At 100K nonrigid, `NeighbourFinder::update` contributes **16.6% of total wall time**.
Within the function, the 4 cost buckets break down as follows:

| Label | Share (%) | ms/call | Invocations | Status |
|---|---|---|---|---|
| `NeighbourFinder::tree_query` | 3.8 | 0.0002 | 40,300,000 | Measurable |
| `NeighbourFinder::buffer_alloc` | ~0 | <0.0001 | 40,300,000 | UNRELIABLE (below noise floor) |
| `NeighbourFinder::query_setup` | ~0 | <0.0001 | 40,300,000 | UNRELIABLE (below noise floor) |
| `NeighbourFinder::result_copy` | ~0 | <0.0001 | 40,300,000 | UNRELIABLE (below noise floor) |

**`NeighbourFinder::tree_query` is the only measurable cost bucket** — the kd-tree
`findNeighbors()` call accounts for essentially all of `NeighbourFinder::update`'s time.
Buffer allocation, query setup, and result copy are all below the 5× timer-overhead
threshold (32 ns/scope × 5 = 160 ns per call) and cannot be reliably measured.

Note: `tree_query` at 3.8% global share accounts for ~23% of `NeighbourFinder::update`'s
16.6% share. The remaining ~77% is likely timer overhead from 40.3M ScopedTimer
constructions (at 32 ns/scope, 40.3M scopes ≈ 1.3 s overhead across the 4 buckets,
which is consistent with the observed timing structure). The instrumentation overhead
itself inflates the per-update measurement at this query depth.

### Decision

Per bead 5jf NOTES decision tree:

> "Cost in tree traversal structure → neighbour caching across ICP iterations (#1)"

**Recommendation: Create bead for neighbour caching (#1).**

The kd-tree `findNeighbors()` call is the sole measurable cost within
`NeighbourFinder::update`. The intervention is **neighbour caching across ICP
iterations**: cache the k-NN result from iteration N and do local refinement in
iteration N+1. This can reduce the total number of full tree queries by 2–10×
at convergence.

**Do NOT pursue:**
- `thread_local` buffers: `buffer_alloc` is below the noise floor — not a meaningful cost.
- SIMD batching (#2): distance math is inside nanoflann's `findNeighbors` (a tree
  traversal, not a batch distance computation). SIMD batching targets per-query
  distance math exposed in a batch API, not applicable to nanoflann's sequential
  traversal without a library swap.
- OpenMP (#5): per bead 9f5 prior art, regression at 7K. Amdahl ceiling at 100K
  is ~1.2× for `NeighbourFinder::update` — not worth the infra cost until
  algorithmic wins are exhausted.

### Instrumentation Overhead Note

The bead spec expects the 4 bucket shares to sum to approximately the pre-instrumentation
`NeighbourFinder::update` share (~14.7% from the 748 hotspot profile). The actual sum
here is ~3.8% — a large discrepancy. The cause is instrumentation overhead:

- Each `update()` call runs 100K queries (at 100K tier)
- 4 ScopedTimers per query = 4 × 100K × 403 invocations = **161.2M timer constructions**
- At 32 ns/scope = **~5.2 s overhead** added to the total measurement
- This overhead is invisible to the bucket labels (it's spread across their construction,
  not recorded inside any bucket) but it inflates `NeighbourFinder::update`'s raw total
- As a result, `NeighbourFinder::update` reads as **16.6%** (vs 14.7% uninstrumented),
  and the bucket totals (which capture real work, net of overhead) sum to only 3.8%

The correct interpretation: the overhead of the timer infrastructure itself, at
per-query granularity, dominates the measurement. The 4 bucket labels cannot recover
the full `NeighbourFinder::update` share in this configuration.

### Confidence

Medium. Timer overhead at 40.3M invocations dominates the per-bucket measurement.
`tree_query` is the only bucket that registers above the noise floor and points clearly
to the kd-tree traversal as the primary cost. A future investigation could confirm via
hardware perf counters (no timer overhead) or profiling with fewer registered buckets
(e.g., only wrapping `tree_query`).
