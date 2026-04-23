# MeshMonk Hotspot Profile — 2026-04-23 (`eb81379`)

## Methodology

- **CPU**: aarch64 (14 logical processors)
- **OS**: Linux 6.17.8-orbstack-00308-g8f9c941121b1
- **Compiler**: cc (Debian 12.2.0-14+deb12u1) 12.2.0
- **CMAKE_BUILD_TYPE**: Release (assumed; use `cmake --build . --config Release`)
- **LTO**: unknown
- **Eigen**: vendored (vendor/eigen-3.4.0)
- **nanoflann**: vendored (vendor/nanoflann.hpp)
- **OpenMesh**: vendored (vendor/OpenMesh-11.0.0)
- **Repo SHA**: `eb813792b6b64af5bdd3a288ef0415357d8427d6`
- **OMP_NUM_THREADS**: 1
- **EIGEN_DONT_PARALLELIZE**: 1
- **Run command**: `python -m profiling.run_profile --tiers 1k,10k,100k --modes rigid,nonrigid,pyramid --runs 5 --warmup 1 --seed 42`
- **Seed**: 42
- **ScopedTimer overhead**: 34 ns/scope (measured via profiling_calibrate(1_000_000))
- **Repetitions per (tier, mode)**: warmup=1, runs=5, median selected

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


### Tier: 1K (~1,000 vertices) — rigid

| Label | Share (%) | ms/call | Invocations | Amdahl ceiling | Fork-join cost (us) | Notes |
|---|---|---|---|---|---|---|
| `RigidRegistration::update` | 35.7 | 70.6380 | 1 | 1.6x | 15 |  |
| `SymmetricCorrespondenceFilter::update` | 26.1 | 0.6374 | 81 | 1.4x | 1215 |  |
| `CorrespondenceFilter::update` | 19.5 | 0.2383 | 162 | 1.2x | 2430 |  |
| `NeighbourFinder::update` | 16.4 | 0.2006 | 162 | 1.2x | 2430 |  |
| `InlierDetector::update` | 2.4 | 0.0580 | 81 | 1.0x | 1215 |  |

*UNRELIABLE = ms/call < 5× timer overhead (34 ns/scope = 0.0000 ms). Excluded from recommendations.*


### Tier: 1K (~1,000 vertices) — nonrigid

| Label | Share (%) | ms/call | Invocations | Amdahl ceiling | Fork-join cost (us) | Notes |
|---|---|---|---|---|---|---|
| `NonrigidRegistration::update` | 41.3 | 284.0790 | 1 | 1.7x | 15 |  |
| `SymmetricCorrespondenceFilter::update` | 16.8 | 0.5734 | 201 | 1.2x | 3015 |  |
| `CorrespondenceFilter::update` | 12.0 | 0.2051 | 402 | 1.1x | 6030 |  |
| `NeighbourFinder::update` | 9.9 | 0.1694 | 403 | 1.1x | 6045 |  |
| `ViscoElasticTransformer::_update_elastically` | 8.5 | 0.2913 | 200 | 1.1x | 3000 |  |
| `ViscoElasticTransformer::_update_viscously` | 8.3 | 0.2851 | 200 | 1.1x | 3000 |  |
| `InlierDetector::update` | 1.7 | 0.0586 | 201 | 1.0x | 3015 |  |
| `ViscoElasticTransformer::_apply_transformation` | 1.1 | 0.0364 | 200 | 1.0x | 3000 |  |
| `ViscoElasticTransformer::_update_outlier_transformation` | 0.4 | 0.0133 | 200 | 1.0x | 3000 |  |
| `ViscoElasticTransformer::_update_smoothing_weights` | 0.0 | 0.0190 | 1 | 1.0x | 15 |  |

*UNRELIABLE = ms/call < 5× timer overhead (34 ns/scope = 0.0000 ms). Excluded from recommendations.*


### Tier: 1K (~1,000 vertices) — pyramid

| Label | Share (%) | ms/call | Invocations | Amdahl ceiling | Fork-join cost (us) | Notes |
|---|---|---|---|---|---|---|
| `PyramidNonrigidRegistration::update` | 23.5 | 82.3330 | 1 | 1.3x | 15 |  |
| `NonrigidRegistration::update` | 21.3 | 24.9023 | 3 | 1.3x | 45 |  |
| `SymmetricCorrespondenceFilter::update` | 10.7 | 0.4112 | 91 | 1.1x | 1365 |  |
| `NonrigidRegistration::update/layer2` | 8.6 | 30.3180 | 1 | 1.1x | 15 |  |
| `CorrespondenceFilter::update` | 7.6 | 0.1457 | 182 | 1.1x | 2730 |  |
| `NonrigidRegistration::update/layer1` | 7.5 | 26.2260 | 1 | 1.1x | 15 |  |
| `NonrigidRegistration::update/layer0` | 7.3 | 25.6750 | 1 | 1.1x | 15 |  |
| `NeighbourFinder::update` | 6.5 | 0.1207 | 188 | 1.1x | 2820 |  |
| `ViscoElasticTransformer::_update_elastically` | 2.4 | 0.0935 | 90 | 1.0x | 1350 |  |
| `ViscoElasticTransformer::_update_viscously` | 2.3 | 0.0914 | 90 | 1.0x | 1350 |  |
| `InlierDetector::update` | 1.1 | 0.0439 | 91 | 1.0x | 1365 |  |
| `ViscoElasticTransformer::_apply_transformation` | 0.7 | 0.0266 | 90 | 1.0x | 1350 |  |
| `ViscoElasticTransformer::_update_outlier_transformation` | 0.5 | 0.0178 | 90 | 1.0x | 1350 |  |
| `ViscoElasticTransformer::_update_smoothing_weights` | 0.0 | 0.0180 | 3 | 1.0x | 45 |  |

  > **SANITY WARNING**: bare `NonrigidRegistration::update` total (74707 us) diverges from layer sum (82219 us) by 9.1% (>5% tolerance).
*UNRELIABLE = ms/call < 5× timer overhead (34 ns/scope = 0.0000 ms). Excluded from recommendations.*


### Tier: 10K (~10,000 vertices) — rigid

| Label | Share (%) | ms/call | Invocations | Amdahl ceiling | Fork-join cost (us) | Notes |
|---|---|---|---|---|---|---|
| `RigidRegistration::update` | 35.2 | 1390.2050 | 1 | 1.5x | 15 |  |
| `SymmetricCorrespondenceFilter::update` | 26.8 | 13.0694 | 81 | 1.4x | 1215 |  |
| `CorrespondenceFilter::update` | 19.2 | 4.6878 | 162 | 1.2x | 2430 |  |
| `NeighbourFinder::update` | 17.6 | 4.2862 | 162 | 1.2x | 2430 |  |
| `InlierDetector::update` | 1.2 | 0.5709 | 81 | 1.0x | 1215 |  |

*UNRELIABLE = ms/call < 5× timer overhead (34 ns/scope = 0.0000 ms). Excluded from recommendations.*


### Tier: 10K (~10,000 vertices) — nonrigid

| Label | Share (%) | ms/call | Invocations | Amdahl ceiling | Fork-join cost (us) | Notes |
|---|---|---|---|---|---|---|
| `NonrigidRegistration::update` | 39.8 | 4441.4930 | 1 | 1.7x | 15 |  |
| `SymmetricCorrespondenceFilter::update` | 20.9 | 11.6037 | 201 | 1.3x | 3015 |  |
| `CorrespondenceFilter::update` | 14.0 | 3.8956 | 402 | 1.2x | 6030 |  |
| `NeighbourFinder::update` | 12.6 | 3.4854 | 403 | 1.1x | 6045 |  |
| `ViscoElasticTransformer::_update_elastically` | 5.3 | 2.9452 | 200 | 1.1x | 3000 |  |
| `ViscoElasticTransformer::_update_viscously` | 5.2 | 2.9080 | 200 | 1.1x | 3000 |  |
| `InlierDetector::update` | 1.0 | 0.5704 | 201 | 1.0x | 3015 |  |
| `ViscoElasticTransformer::_apply_transformation` | 0.8 | 0.4219 | 200 | 1.0x | 3000 |  |
| `ViscoElasticTransformer::_update_outlier_transformation` | 0.3 | 0.1715 | 200 | 1.0x | 3000 |  |
| `ViscoElasticTransformer::_update_smoothing_weights` | 0.0 | 0.2370 | 1 | 1.0x | 15 |  |

*UNRELIABLE = ms/call < 5× timer overhead (34 ns/scope = 0.0000 ms). Excluded from recommendations.*


### Tier: 10K (~10,000 vertices) — pyramid

| Label | Share (%) | ms/call | Invocations | Amdahl ceiling | Fork-join cost (us) | Notes |
|---|---|---|---|---|---|---|
| `PyramidNonrigidRegistration::update` | 23.0 | 1244.1080 | 1 | 1.3x | 15 |  |
| `NonrigidRegistration::update` | 21.4 | 385.8496 | 3 | 1.3x | 45 |  |
| `SymmetricCorrespondenceFilter::update` | 12.4 | 7.3743 | 91 | 1.1x | 1365 |  |
| `NonrigidRegistration::update/layer2` | 9.2 | 498.5570 | 1 | 1.1x | 15 |  |
| `CorrespondenceFilter::update` | 8.1 | 2.4049 | 182 | 1.1x | 2730 |  |
| `NonrigidRegistration::update/layer1` | 7.4 | 402.4080 | 1 | 1.1x | 15 |  |
| `NeighbourFinder::update` | 7.4 | 2.1205 | 188 | 1.1x | 2820 |  |
| `NonrigidRegistration::update/layer0` | 6.3 | 339.6350 | 1 | 1.1x | 15 |  |
| `ViscoElasticTransformer::_update_elastically` | 1.6 | 0.9609 | 90 | 1.0x | 1350 |  |
| `ViscoElasticTransformer::_update_viscously` | 1.6 | 0.9499 | 90 | 1.0x | 1350 |  |
| `InlierDetector::update` | 0.7 | 0.4324 | 91 | 1.0x | 1365 |  |
| `ViscoElasticTransformer::_apply_transformation` | 0.5 | 0.3196 | 90 | 1.0x | 1350 |  |
| `ViscoElasticTransformer::_update_outlier_transformation` | 0.3 | 0.1558 | 90 | 1.0x | 1350 |  |
| `ViscoElasticTransformer::_update_smoothing_weights` | 0.0 | 0.1743 | 3 | 1.0x | 45 |  |

  > **SANITY WARNING**: bare `NonrigidRegistration::update` total (1157549 us) diverges from layer sum (1240600 us) by 6.7% (>5% tolerance).
*UNRELIABLE = ms/call < 5× timer overhead (34 ns/scope = 0.0000 ms). Excluded from recommendations.*


### Tier: 100K (~100,000 vertices) — rigid

| Label | Share (%) | ms/call | Invocations | Amdahl ceiling | Fork-join cost (us) | Notes |
|---|---|---|---|---|---|---|
| `RigidRegistration::update` | 32.8 | 19321.4900 | 1 | 1.5x | 15 |  |
| `SymmetricCorrespondenceFilter::update` | 26.6 | 193.3635 | 81 | 1.4x | 1215 |  |
| `CorrespondenceFilter::update` | 20.5 | 74.5132 | 162 | 1.3x | 2430 |  |
| `NeighbourFinder::update` | 19.2 | 69.9090 | 162 | 1.2x | 2430 |  |
| `InlierDetector::update` | 0.8 | 5.8911 | 81 | 1.0x | 1215 |  |

*UNRELIABLE = ms/call < 5× timer overhead (34 ns/scope = 0.0000 ms). Excluded from recommendations.*


### Tier: 100K (~100,000 vertices) — nonrigid

| Label | Share (%) | ms/call | Invocations | Amdahl ceiling | Fork-join cost (us) | Notes |
|---|---|---|---|---|---|---|
| `NonrigidRegistration::update` | 37.2 | 56541.1430 | 1 | 1.6x | 15 |  |
| `SymmetricCorrespondenceFilter::update` | 22.2 | 167.8699 | 201 | 1.3x | 3015 |  |
| `CorrespondenceFilter::update` | 16.1 | 60.7200 | 402 | 1.2x | 6030 |  |
| `NeighbourFinder::update` | 14.8 | 55.7294 | 403 | 1.2x | 6045 |  |
| `ViscoElasticTransformer::_update_elastically` | 4.1 | 31.3213 | 200 | 1.0x | 3000 |  |
| `ViscoElasticTransformer::_update_viscously` | 4.1 | 30.8796 | 200 | 1.0x | 3000 |  |
| `InlierDetector::update` | 0.8 | 5.9024 | 201 | 1.0x | 3015 |  |
| `ViscoElasticTransformer::_apply_transformation` | 0.6 | 4.5452 | 200 | 1.0x | 3000 |  |
| `ViscoElasticTransformer::_update_outlier_transformation` | 0.2 | 1.7174 | 200 | 1.0x | 3000 |  |
| `ViscoElasticTransformer::_update_smoothing_weights` | 0.0 | 1.9800 | 1 | 1.0x | 15 |  |

*UNRELIABLE = ms/call < 5× timer overhead (34 ns/scope = 0.0000 ms). Excluded from recommendations.*


### Tier: 100K (~100,000 vertices) — pyramid

| Label | Share (%) | ms/call | Invocations | Amdahl ceiling | Fork-join cost (us) | Notes |
|---|---|---|---|---|---|---|
| `PyramidNonrigidRegistration::update` | 22.1 | 16087.2420 | 1 | 1.3x | 15 |  |
| `NonrigidRegistration::update` | 20.7 | 5024.2636 | 3 | 1.3x | 45 |  |
| `SymmetricCorrespondenceFilter::update` | 13.4 | 107.3969 | 91 | 1.2x | 1365 |  |
| `CorrespondenceFilter::update` | 9.4 | 37.6529 | 182 | 1.1x | 2730 |  |
| `NeighbourFinder::update` | 8.7 | 33.8844 | 188 | 1.1x | 2820 |  |
| `NonrigidRegistration::update/layer2` | 8.7 | 6345.2060 | 1 | 1.1x | 15 |  |
| `NonrigidRegistration::update/layer1` | 7.3 | 5326.1720 | 1 | 1.1x | 15 |  |
| `NonrigidRegistration::update/layer0` | 6.0 | 4395.1280 | 1 | 1.1x | 15 |  |
| `ViscoElasticTransformer::_update_viscously` | 1.3 | 10.3733 | 90 | 1.0x | 1350 |  |
| `ViscoElasticTransformer::_update_elastically` | 1.3 | 10.2597 | 90 | 1.0x | 1350 |  |
| `InlierDetector::update` | 0.6 | 4.4097 | 91 | 1.0x | 1365 |  |
| `ViscoElasticTransformer::_apply_transformation` | 0.4 | 3.2030 | 90 | 1.0x | 1350 |  |
| `ViscoElasticTransformer::_update_outlier_transformation` | 0.2 | 1.4515 | 90 | 1.0x | 1350 |  |
| `ViscoElasticTransformer::_update_smoothing_weights` | 0.0 | 1.8520 | 3 | 1.0x | 45 |  |

  > **SANITY WARNING**: bare `NonrigidRegistration::update` total (15072791 us) diverges from layer sum (16066506 us) by 6.2% (>5% tolerance).
*UNRELIABLE = ms/call < 5× timer overhead (34 ns/scope = 0.0000 ms). Excluded from recommendations.*

## Scaling Exponents (log-log fit across 1K / 10K / 100K)

| Label | Mode | k (exponent) | Complexity bucket | Note |
|---|---|---|---|---|
| `InlierDetector::update` | rigid | 1.00 | ≈O(n^1.0) |  |
| `CorrespondenceFilter::update` | rigid | 1.25 | ≈O(n^1.16) |  |
| `SymmetricCorrespondenceFilter::update` | rigid | 1.24 | ≈O(n^1.16) |  |
| `NeighbourFinder::update` | rigid | 1.27 | ≈O(n^1.16) |  |
| `RigidRegistration::update` | rigid | 1.22 | ≈O(n^1.16) |  |
| `InlierDetector::update` | nonrigid | 1.00 | ≈O(n^1.0) |  |
| `CorrespondenceFilter::update` | nonrigid | 1.24 | ≈O(n^1.16) |  |
| `ViscoElasticTransformer::_update_viscously` | nonrigid | 1.02 | ≈O(n^1.0) |  |
| `ViscoElasticTransformer::_update_outlier_transformation` | nonrigid | 1.06 | ≈O(n^1.0) |  |
| `ViscoElasticTransformer::_apply_transformation` | nonrigid | 1.05 | ≈O(n^1.0) |  |
| `SymmetricCorrespondenceFilter::update` | nonrigid | 1.23 | ≈O(n^1.16) |  |
| `NeighbourFinder::update` | nonrigid | 1.26 | ≈O(n^1.16) |  |
| `ViscoElasticTransformer::_update_elastically` | nonrigid | 1.02 | ≈O(n^1.0) |  |
| `NonrigidRegistration::update` | nonrigid | 1.15 | ≈O(n^1.0) |  |
| `ViscoElasticTransformer::_update_smoothing_weights` | nonrigid | 1.01 | ≈O(n^1.0) |  |
| `InlierDetector::update` | pyramid | 1.00 | ≈O(n^1.0) |  |
| `NonrigidRegistration::update/layer1` | pyramid | 1.15 | ≈O(n^1.0) |  |
| `PyramidNonrigidRegistration::update` | pyramid | 1.15 | ≈O(n^1.0) |  |
| `ViscoElasticTransformer::_update_smoothing_weights` | pyramid | 1.01 | ≈O(n^1.0) |  |
| `NonrigidRegistration::update/layer2` | pyramid | 1.16 | ≈O(n^1.0) |  |
| `ViscoElasticTransformer::_apply_transformation` | pyramid | 1.04 | ≈O(n^1.0) |  |
| `ViscoElasticTransformer::_update_outlier_transformation` | pyramid | 0.96 | ≈O(n^1.0) |  |
| `CorrespondenceFilter::update` | pyramid | 1.21 | ≈O(n^1.16) |  |
| `NonrigidRegistration::update/layer0` | pyramid | 1.12 | ≈O(n^1.0) |  |
| `SymmetricCorrespondenceFilter::update` | pyramid | 1.21 | ≈O(n^1.16) |  |
| `NeighbourFinder::update` | pyramid | 1.22 | ≈O(n^1.16) |  |
| `ViscoElasticTransformer::_update_elastically` | pyramid | 1.02 | ≈O(n^1.0) |  |
| `NonrigidRegistration::update` | pyramid | 1.15 | ≈O(n^1.0) |  |
| `ViscoElasticTransformer::_update_viscously` | pyramid | 1.03 | ≈O(n^1.0) |  |

> **Caveat:** Three-point fits (1K/10K/100K) have zero residual degrees of freedom.
> The buckets k≈1.0 / k≈1.16 / k≈2.0 cannot be reliably discriminated without a
> fourth tier or per-tier variance. Treat exponents as indicative only.

## Recommendations

### Candidates for optimization

Labels meeting all of: share > 5%, ms/call > 0.03 ms, NOT UNRELIABLE.

- **`NonrigidRegistration::update`** (nonrigid, 1K tier): share=41.3%, Amdahl ceiling = 1/(1−0.413) = 1.7×
  - Scaling exponent k: 1.15
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 15 µs (0.0% of label wall time)

- **`RigidRegistration::update`** (rigid, 1K tier): share=35.7%, Amdahl ceiling = 1/(1−0.357) = 1.6×
  - Scaling exponent k: 1.22
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 15 µs (0.0% of label wall time)

- **`SymmetricCorrespondenceFilter::update`** (rigid, 10K tier): share=26.8%, Amdahl ceiling = 1/(1−0.268) = 1.4×
  - Scaling exponent k: 1.24
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 1215 µs (0.1% of label wall time)

- **`PyramidNonrigidRegistration::update`** (pyramid, 1K tier): share=23.5%, Amdahl ceiling = 1/(1−0.235) = 1.3×
  - Scaling exponent k: 1.15
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 15 µs (0.0% of label wall time)

- **`SymmetricCorrespondenceFilter::update`** (nonrigid, 100K tier): share=22.2%, Amdahl ceiling = 1/(1−0.222) = 1.3×
  - Scaling exponent k: 1.23
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 3015 µs (0.0% of label wall time)

- **`NonrigidRegistration::update`** (pyramid, 10K tier): share=21.4%, Amdahl ceiling = 1/(1−0.214) = 1.3×
  - Scaling exponent k: 1.15
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 45 µs (0.0% of label wall time)

- **`CorrespondenceFilter::update`** (rigid, 100K tier): share=20.5%, Amdahl ceiling = 1/(1−0.205) = 1.3×
  - Scaling exponent k: 1.25
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 2430 µs (0.0% of label wall time)

- **`NeighbourFinder::update`** (rigid, 100K tier): share=19.2%, Amdahl ceiling = 1/(1−0.192) = 1.2×
  - Scaling exponent k: 1.27
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 2430 µs (0.0% of label wall time)

- **`CorrespondenceFilter::update`** (nonrigid, 100K tier): share=16.1%, Amdahl ceiling = 1/(1−0.161) = 1.2×
  - Scaling exponent k: 1.24
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 6030 µs (0.0% of label wall time)

- **`NeighbourFinder::update`** (nonrigid, 100K tier): share=14.8%, Amdahl ceiling = 1/(1−0.148) = 1.2×
  - Scaling exponent k: 1.26
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 6045 µs (0.0% of label wall time)

- **`SymmetricCorrespondenceFilter::update`** (pyramid, 100K tier): share=13.4%, Amdahl ceiling = 1/(1−0.134) = 1.2×
  - Scaling exponent k: 1.21
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 1365 µs (0.0% of label wall time)

- **`CorrespondenceFilter::update`** (pyramid, 100K tier): share=9.4%, Amdahl ceiling = 1/(1−0.094) = 1.1×
  - Scaling exponent k: 1.21
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 2730 µs (0.0% of label wall time)

- **`NonrigidRegistration::update/layer2`** (pyramid, 10K tier): share=9.2%, Amdahl ceiling = 1/(1−0.092) = 1.1×
  - Scaling exponent k: 1.16
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 15 µs (0.0% of label wall time)

- **`NeighbourFinder::update`** (pyramid, 100K tier): share=8.7%, Amdahl ceiling = 1/(1−0.087) = 1.1×
  - Scaling exponent k: 1.22
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 2820 µs (0.0% of label wall time)

- **`ViscoElasticTransformer::_update_elastically`** (nonrigid, 1K tier): share=8.5%, Amdahl ceiling = 1/(1−0.085) = 1.1×
  - Scaling exponent k: 1.02
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 3000 µs (5.1% of label wall time)

- **`ViscoElasticTransformer::_update_viscously`** (nonrigid, 1K tier): share=8.3%, Amdahl ceiling = 1/(1−0.083) = 1.1×
  - Scaling exponent k: 1.02
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 3000 µs (5.3% of label wall time)

- **`NonrigidRegistration::update/layer1`** (pyramid, 1K tier): share=7.5%, Amdahl ceiling = 1/(1−0.075) = 1.1×
  - Scaling exponent k: 1.15
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 15 µs (0.1% of label wall time)

- **`NonrigidRegistration::update/layer0`** (pyramid, 1K tier): share=7.3%, Amdahl ceiling = 1/(1−0.073) = 1.1×
  - Scaling exponent k: 1.12
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 15 µs (0.1% of label wall time)

### Below granularity threshold

Labels with share ≤ 5% or ms/call ≤ 0.03 ms. OpenMP not viable due to
fork-join overhead exceeding potential speedup (reference: beads 9f5, bdt).

`InlierDetector::update`, `ViscoElasticTransformer::_apply_transformation`, `ViscoElasticTransformer::_update_elastically`, `ViscoElasticTransformer::_update_outlier_transformation`, `ViscoElasticTransformer::_update_smoothing_weights`, `ViscoElasticTransformer::_update_viscously`

### Excluded (UNRELIABLE)

Labels with ms/call < 5× timer overhead (34 ns/scope = 0.0000 ms).
Corrected values unreliable — cannot meaningfully distinguish measurement noise from actual cost.

(none)
