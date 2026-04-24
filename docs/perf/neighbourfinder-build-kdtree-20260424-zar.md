# MeshMonk Hotspot Profile — 2026-04-24 (`9dd5bf8`)

## Methodology

- **CPU**: arm (14 logical processors)
- **OS**: Darwin 24.6.0
- **Compiler**: Apple clang version 17.0.0 (clang-1700.0.13.5)
- **CMAKE_BUILD_TYPE**: Release (assumed; use `cmake --build . --config Release`)
- **LTO**: unknown
- **Eigen**: vendored (vendor/eigen-3.4.0)
- **nanoflann**: vendored (vendor/nanoflann.hpp)
- **OpenMesh**: vendored (vendor/OpenMesh-11.0.0)
- **Repo SHA**: `9dd5bf833c7eebcbb5e7daa6cb5ea048a9daffaa`
- **OMP_NUM_THREADS**: 1
- **EIGEN_DONT_PARALLELIZE**: 1
- **Run command**: `python -m profiling.run_profile --tiers 100k --modes nonrigid --runs 3 --warmup 1 --seed 42`
- **Seed**: 42
- **ScopedTimer overhead**: 35 ns/scope (median of 3 × profiling_calibrate(1_000_000); individual: 35, 35, 36 ns/scope)
- **Repetitions per (tier, mode)**: warmup=1, runs=3, median selected

### Mesh SHA-256

| File | SHA-256 |
|---|---|
| Template_1K.obj | `0abb45cb88d1f204577826e4f614f008fc288d4424215f1d08d15d923f6849c5` |
| Template_10K.obj | `bd8ceb706e5e0ecf19eff74f40bde6459f3a7061a5c777b7fe6989cbf1fef087` |
| Template_100K.obj | `d2353cae5fbdf6216dfc918ac3f28932e8a766ad38c94ce739b1a67f9e8157f4` |
| DemoFace_1K.obj | `7f8ec6a5b17671750a135717cb223c533f5a858052ffb9cc007c22ed76494da6` |
| DemoFace_10K.obj | `460a5109ce8e8b00a704796b5e05859058662dd255538fe2595f2cf4852d9e4b` |
| DemoFace_100K.obj | `26610e15f32b70422adde365f46fcb2b5e3296da21afca3091b21c636aa5922e` |

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
| `NonrigidRegistration::update` | 33.1 | 62654.4910 | 1 | 1.5x | 15 |  |
| `SymmetricCorrespondenceFilter::update` | 20.0 | 188.3079 | 201 | 1.3x | 3015 |  |
| `CorrespondenceFilter::update` | 13.7 | 64.6120 | 402 | 1.2x | 6030 |  |
| `NeighbourFinder::update` | 12.6 | 59.2748 | 403 | 1.1x | 6045 |  |
| `NeighbourFinder::build_kdtree` | 10.7 | 25.0943 | 807 | 1.1x | 12105 |  |
| `ViscoElasticTransformer::_update_viscously` | 3.3 | 31.2840 | 200 | 1.0x | 3000 |  |
| `ViscoElasticTransformer::_update_elastically` | 3.3 | 31.2572 | 200 | 1.0x | 3000 |  |
| `NeighbourFinder::tree_query` | 2.0 | 0.0001 | 40300000 | 1.0x | 604500000 | * UNRELIABLE |
| `InlierDetector::update` | 0.6 | 5.3033 | 201 | 1.0x | 3015 |  |
| `ViscoElasticTransformer::_apply_transformation` | 0.4 | 3.7525 | 200 | 1.0x | 3000 |  |
| `ViscoElasticTransformer::_update_outlier_transformation` | 0.2 | 1.7144 | 200 | 1.0x | 3000 |  |
| `ViscoElasticTransformer::_update_smoothing_weights` | 0.0 | 2.6380 | 1 | 1.0x | 15 |  |

*UNRELIABLE = ms/call < 5× timer overhead (35 ns/scope = 0.000035 ms). Excluded from recommendations.*

## Scaling Exponents (log-log fit across 1K / 10K / 100K)

| Label | Mode | k (exponent) | Complexity bucket | Note |
|---|---|---|---|---|
| `ViscoElasticTransformer::_update_smoothing_weights` | nonrigid | NaN | n/a (insufficient tiers) | insufficient tiers |
| `NeighbourFinder::tree_query` | nonrigid | NaN | n/a (insufficient tiers) | insufficient tiers |
| `NeighbourFinder::update` | nonrigid | NaN | n/a (insufficient tiers) | insufficient tiers |
| `ViscoElasticTransformer::_update_viscously` | nonrigid | NaN | n/a (insufficient tiers) | insufficient tiers |
| `SymmetricCorrespondenceFilter::update` | nonrigid | NaN | n/a (insufficient tiers) | insufficient tiers |
| `NeighbourFinder::build_kdtree` | nonrigid | NaN | n/a (insufficient tiers) | insufficient tiers |
| `ViscoElasticTransformer::_update_outlier_transformation` | nonrigid | NaN | n/a (insufficient tiers) | insufficient tiers |
| `ViscoElasticTransformer::_update_elastically` | nonrigid | NaN | n/a (insufficient tiers) | insufficient tiers |
| `CorrespondenceFilter::update` | nonrigid | NaN | n/a (insufficient tiers) | insufficient tiers |
| `NonrigidRegistration::update` | nonrigid | NaN | n/a (insufficient tiers) | insufficient tiers |
| `ViscoElasticTransformer::_apply_transformation` | nonrigid | NaN | n/a (insufficient tiers) | insufficient tiers |
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

- **`NonrigidRegistration::update`** (nonrigid, 100K tier): share=33.1%, Amdahl ceiling = 1/(1−0.331) = 1.5×
  - Scaling exponent k: NaN
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 15 µs (0.0% of label wall time)

- **`SymmetricCorrespondenceFilter::update`** (nonrigid, 100K tier): share=20.0%, Amdahl ceiling = 1/(1−0.200) = 1.3×
  - Scaling exponent k: NaN
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 3015 µs (0.0% of label wall time)

- **`CorrespondenceFilter::update`** (nonrigid, 100K tier): share=13.7%, Amdahl ceiling = 1/(1−0.137) = 1.2×
  - Scaling exponent k: NaN
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 6030 µs (0.0% of label wall time)

- **`NeighbourFinder::update`** (nonrigid, 100K tier): share=12.6%, Amdahl ceiling = 1/(1−0.126) = 1.1×
  - Scaling exponent k: NaN
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 6045 µs (0.0% of label wall time)

- **`NeighbourFinder::build_kdtree`** (nonrigid, 100K tier): share=10.7%, Amdahl ceiling = 1/(1−0.107) = 1.1×
  - Scaling exponent k: NaN
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 12105 µs (0.1% of label wall time)

### Below granularity threshold

Labels with share ≤ 5% or ms/call ≤ 0.03 ms. OpenMP not viable due to
fork-join overhead exceeding potential speedup (reference: beads 9f5, bdt).

`InlierDetector::update`, `ViscoElasticTransformer::_apply_transformation`, `ViscoElasticTransformer::_update_elastically`, `ViscoElasticTransformer::_update_outlier_transformation`, `ViscoElasticTransformer::_update_smoothing_weights`, `ViscoElasticTransformer::_update_viscously`

### Excluded (UNRELIABLE)

Labels with ms/call < 5× timer overhead (35 ns/scope = 0.000035 ms).
Corrected values unreliable — cannot meaningfully distinguish measurement noise from actual cost.

| Label | Mode | Tier | ms/call |
|---|---|---|---|
| `NeighbourFinder::tree_query` | nonrigid | 100K | 0.0001 |
