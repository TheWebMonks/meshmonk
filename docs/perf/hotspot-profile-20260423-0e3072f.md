# MeshMonk Hotspot Profile — 2026-04-23 (`0e3072f`)

## Methodology

- **CPU**: aarch64 (14 logical processors)
- **OS**: Linux 6.17.8-orbstack-00308-g8f9c941121b1
- **Compiler**: cc (Debian 12.2.0-14+deb12u1) 12.2.0
- **CMAKE_BUILD_TYPE**: Release (assumed; use `cmake --build . --config Release`)
- **LTO**: unknown
- **Eigen**: vendored (vendor/eigen-3.4.0)
- **nanoflann**: vendored (vendor/nanoflann.hpp)
- **OpenMesh**: vendored (vendor/OpenMesh-11.0.0)
- **Repo SHA**: `0e3072fd60d4e5913cf8effe7ddf0997e33caed0`
- **OMP_NUM_THREADS**: 1
- **EIGEN_DONT_PARALLELIZE**: 1
- **Run command**: `python -m profiling.run_profile --tiers 1k,10k,100k --modes rigid,nonrigid,pyramid --runs 5 --warmup 1 --seed 42`
- **Seed**: 42
- **ScopedTimer overhead**: 25 ns/scope (median of 3 × profiling_calibrate(1_000_000); individual: 25, 26, 25 ns/scope)
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
| `RigidRegistration::update` | 35.5 | 71.3220 | 1 | 1.6x | 15 |  |
| `SymmetricCorrespondenceFilter::update` | 26.0 | 0.6458 | 81 | 1.4x | 1215 |  |
| `CorrespondenceFilter::update` | 19.6 | 0.2430 | 162 | 1.2x | 2430 |  |
| `NeighbourFinder::update` | 16.5 | 0.2052 | 162 | 1.2x | 2430 |  |
| `InlierDetector::update` | 2.3 | 0.0580 | 81 | 1.0x | 1215 |  |

*UNRELIABLE = ms/call < 5× timer overhead (25 ns/scope = 0.000025 ms). Excluded from recommendations.*


### Tier: 1K (~1,000 vertices) — nonrigid

| Label | Share (%) | ms/call | Invocations | Amdahl ceiling | Fork-join cost (us) | Notes |
|---|---|---|---|---|---|---|
| `NonrigidRegistration::update` | 41.3 | 284.9950 | 1 | 1.7x | 15 |  |
| `SymmetricCorrespondenceFilter::update` | 16.9 | 0.5801 | 201 | 1.2x | 3015 |  |
| `CorrespondenceFilter::update` | 12.1 | 0.2079 | 402 | 1.1x | 6030 |  |
| `NeighbourFinder::update` | 10.0 | 0.1721 | 403 | 1.1x | 6045 |  |
| `ViscoElasticTransformer::_update_elastically` | 8.3 | 0.2881 | 200 | 1.1x | 3000 |  |
| `ViscoElasticTransformer::_update_viscously` | 8.2 | 0.2827 | 200 | 1.1x | 3000 |  |
| `InlierDetector::update` | 1.7 | 0.0589 | 201 | 1.0x | 3015 |  |
| `ViscoElasticTransformer::_apply_transformation` | 1.1 | 0.0365 | 200 | 1.0x | 3000 |  |
| `ViscoElasticTransformer::_update_outlier_transformation` | 0.4 | 0.0136 | 200 | 1.0x | 3000 |  |
| `ViscoElasticTransformer::_update_smoothing_weights` | 0.0 | 0.0200 | 1 | 1.0x | 15 |  |

*UNRELIABLE = ms/call < 5× timer overhead (25 ns/scope = 0.000025 ms). Excluded from recommendations.*


### Tier: 1K (~1,000 vertices) — pyramid

| Label | Share (%) | ms/call | Invocations | Amdahl ceiling | Fork-join cost (us) | Notes |
|---|---|---|---|---|---|---|
| `PyramidNonrigidRegistration::update` | 24.0 | 82.4020 | 1 | 1.3x | 15 |  |
| `NonrigidRegistration::update` | 21.8 | 24.9476 | 3 | 1.3x | 45 |  |
| `SymmetricCorrespondenceFilter::update` | 10.9 | 0.4123 | 91 | 1.1x | 1365 |  |
| `NonrigidRegistration::update/layer2` | 8.6 | 29.4770 | 1 | 1.1x | 15 |  |
| `CorrespondenceFilter::update` | 7.8 | 0.1468 | 182 | 1.1x | 2730 |  |
| `NonrigidRegistration::update/layer1` | 6.8 | 23.3550 | 1 | 1.1x | 15 |  |
| `NeighbourFinder::update` | 6.7 | 0.1219 | 188 | 1.1x | 2820 |  |
| `NonrigidRegistration::update/layer0` | 6.4 | 22.0260 | 1 | 1.1x | 15 |  |
| `ViscoElasticTransformer::_update_elastically` | 2.4 | 0.0932 | 90 | 1.0x | 1350 |  |
| `ViscoElasticTransformer::_update_viscously` | 2.4 | 0.0922 | 90 | 1.0x | 1350 |  |
| `InlierDetector::update` | 1.2 | 0.0439 | 91 | 1.0x | 1365 |  |
| `ViscoElasticTransformer::_apply_transformation` | 0.7 | 0.0266 | 90 | 1.0x | 1350 |  |
| `ViscoElasticTransformer::_update_outlier_transformation` | 0.5 | 0.0184 | 90 | 1.0x | 1350 |  |
| `ViscoElasticTransformer::_update_smoothing_weights` | 0.0 | 0.0166 | 3 | 1.0x | 45 |  |

*UNRELIABLE = ms/call < 5× timer overhead (25 ns/scope = 0.000025 ms). Excluded from recommendations.*


### Tier: 10K (~10,000 vertices) — rigid

| Label | Share (%) | ms/call | Invocations | Amdahl ceiling | Fork-join cost (us) | Notes |
|---|---|---|---|---|---|---|
| `RigidRegistration::update` | 35.3 | 1403.1710 | 1 | 1.5x | 15 |  |
| `SymmetricCorrespondenceFilter::update` | 26.8 | 13.1796 | 81 | 1.4x | 1215 |  |
| `CorrespondenceFilter::update` | 19.2 | 4.7069 | 162 | 1.2x | 2430 |  |
| `NeighbourFinder::update` | 17.5 | 4.3032 | 162 | 1.2x | 2430 |  |
| `InlierDetector::update` | 1.2 | 0.5704 | 81 | 1.0x | 1215 |  |

*UNRELIABLE = ms/call < 5× timer overhead (25 ns/scope = 0.000025 ms). Excluded from recommendations.*


### Tier: 10K (~10,000 vertices) — nonrigid

| Label | Share (%) | ms/call | Invocations | Amdahl ceiling | Fork-join cost (us) | Notes |
|---|---|---|---|---|---|---|
| `NonrigidRegistration::update` | 39.9 | 4487.9020 | 1 | 1.7x | 15 |  |
| `SymmetricCorrespondenceFilter::update` | 21.0 | 11.7401 | 201 | 1.3x | 3015 |  |
| `CorrespondenceFilter::update` | 14.0 | 3.9184 | 402 | 1.2x | 6030 |  |
| `NeighbourFinder::update` | 12.6 | 3.4995 | 403 | 1.1x | 6045 |  |
| `ViscoElasticTransformer::_update_elastically` | 5.2 | 2.9367 | 200 | 1.1x | 3000 |  |
| `ViscoElasticTransformer::_update_viscously` | 5.2 | 2.8996 | 200 | 1.1x | 3000 |  |
| `InlierDetector::update` | 1.0 | 0.5784 | 201 | 1.0x | 3015 |  |
| `ViscoElasticTransformer::_apply_transformation` | 0.8 | 0.4249 | 200 | 1.0x | 3000 |  |
| `ViscoElasticTransformer::_update_outlier_transformation` | 0.3 | 0.1759 | 200 | 1.0x | 3000 |  |
| `ViscoElasticTransformer::_update_smoothing_weights` | 0.0 | 0.2550 | 1 | 1.0x | 15 |  |

*UNRELIABLE = ms/call < 5× timer overhead (25 ns/scope = 0.000025 ms). Excluded from recommendations.*


### Tier: 10K (~10,000 vertices) — pyramid

| Label | Share (%) | ms/call | Invocations | Amdahl ceiling | Fork-join cost (us) | Notes |
|---|---|---|---|---|---|---|
| `PyramidNonrigidRegistration::update` | 23.4 | 1254.9200 | 1 | 1.3x | 15 |  |
| `NonrigidRegistration::update` | 21.8 | 389.5946 | 3 | 1.3x | 45 |  |
| `SymmetricCorrespondenceFilter::update` | 12.6 | 7.4359 | 91 | 1.1x | 1365 |  |
| `NonrigidRegistration::update/layer2` | 9.3 | 497.3100 | 1 | 1.1x | 15 |  |
| `CorrespondenceFilter::update` | 8.2 | 2.4144 | 182 | 1.1x | 2730 |  |
| `NeighbourFinder::update` | 7.4 | 2.1239 | 188 | 1.1x | 2820 |  |
| `NonrigidRegistration::update/layer1` | 7.0 | 374.6210 | 1 | 1.1x | 15 |  |
| `NonrigidRegistration::update/layer0` | 5.5 | 296.2170 | 1 | 1.1x | 15 |  |
| `ViscoElasticTransformer::_update_elastically` | 1.6 | 0.9656 | 90 | 1.0x | 1350 |  |
| `ViscoElasticTransformer::_update_viscously` | 1.6 | 0.9597 | 90 | 1.0x | 1350 |  |
| `InlierDetector::update` | 0.7 | 0.4304 | 91 | 1.0x | 1365 |  |
| `ViscoElasticTransformer::_apply_transformation` | 0.5 | 0.3236 | 90 | 1.0x | 1350 |  |
| `ViscoElasticTransformer::_update_outlier_transformation` | 0.3 | 0.1642 | 90 | 1.0x | 1350 |  |
| `ViscoElasticTransformer::_update_smoothing_weights` | 0.0 | 0.1743 | 3 | 1.0x | 45 |  |

*UNRELIABLE = ms/call < 5× timer overhead (25 ns/scope = 0.000025 ms). Excluded from recommendations.*


### Tier: 100K (~100,000 vertices) — rigid

| Label | Share (%) | ms/call | Invocations | Amdahl ceiling | Fork-join cost (us) | Notes |
|---|---|---|---|---|---|---|
| `RigidRegistration::update` | 32.9 | 19581.2500 | 1 | 1.5x | 15 |  |
| `SymmetricCorrespondenceFilter::update` | 26.6 | 195.5057 | 81 | 1.4x | 1215 |  |
| `CorrespondenceFilter::update` | 20.4 | 75.0312 | 162 | 1.3x | 2430 |  |
| `NeighbourFinder::update` | 19.2 | 70.3893 | 162 | 1.2x | 2430 |  |
| `InlierDetector::update` | 0.8 | 5.8998 | 81 | 1.0x | 1215 |  |

*UNRELIABLE = ms/call < 5× timer overhead (25 ns/scope = 0.000025 ms). Excluded from recommendations.*


### Tier: 100K (~100,000 vertices) — nonrigid

| Label | Share (%) | ms/call | Invocations | Amdahl ceiling | Fork-join cost (us) | Notes |
|---|---|---|---|---|---|---|
| `NonrigidRegistration::update` | 37.3 | 57365.1150 | 1 | 1.6x | 15 |  |
| `SymmetricCorrespondenceFilter::update` | 22.1 | 169.3601 | 201 | 1.3x | 3015 |  |
| `CorrespondenceFilter::update` | 16.0 | 61.0805 | 402 | 1.2x | 6030 |  |
| `NeighbourFinder::update` | 14.7 | 56.0491 | 403 | 1.2x | 6045 |  |
| `ViscoElasticTransformer::_update_elastically` | 4.2 | 32.4521 | 200 | 1.0x | 3000 |  |
| `ViscoElasticTransformer::_update_viscously` | 4.1 | 31.8775 | 200 | 1.0x | 3000 |  |
| `InlierDetector::update` | 0.8 | 5.9027 | 201 | 1.0x | 3015 |  |
| `ViscoElasticTransformer::_apply_transformation` | 0.6 | 4.6453 | 200 | 1.0x | 3000 |  |
| `ViscoElasticTransformer::_update_outlier_transformation` | 0.2 | 1.7718 | 200 | 1.0x | 3000 |  |
| `ViscoElasticTransformer::_update_smoothing_weights` | 0.0 | 2.5400 | 1 | 1.0x | 15 |  |

*UNRELIABLE = ms/call < 5× timer overhead (25 ns/scope = 0.000025 ms). Excluded from recommendations.*


### Tier: 100K (~100,000 vertices) — pyramid

| Label | Share (%) | ms/call | Invocations | Amdahl ceiling | Fork-join cost (us) | Notes |
|---|---|---|---|---|---|---|
| `PyramidNonrigidRegistration::update` | 22.5 | 16717.3730 | 1 | 1.3x | 15 |  |
| `NonrigidRegistration::update` | 21.0 | 5200.3913 | 3 | 1.3x | 45 |  |
| `SymmetricCorrespondenceFilter::update` | 13.6 | 110.8636 | 91 | 1.2x | 1365 |  |
| `CorrespondenceFilter::update` | 9.5 | 38.6815 | 182 | 1.1x | 2730 |  |
| `NonrigidRegistration::update/layer2` | 8.8 | 6544.5270 | 1 | 1.1x | 15 |  |
| `NeighbourFinder::update` | 8.8 | 34.7490 | 188 | 1.1x | 2820 |  |
| `NonrigidRegistration::update/layer1` | 6.9 | 5101.9840 | 1 | 1.1x | 15 |  |
| `NonrigidRegistration::update/layer0` | 5.4 | 3984.2920 | 1 | 1.1x | 15 |  |
| `ViscoElasticTransformer::_update_viscously` | 1.3 | 10.6890 | 90 | 1.0x | 1350 |  |
| `ViscoElasticTransformer::_update_elastically` | 1.3 | 10.6465 | 90 | 1.0x | 1350 |  |
| `InlierDetector::update` | 0.5 | 4.4640 | 91 | 1.0x | 1365 |  |
| `ViscoElasticTransformer::_apply_transformation` | 0.4 | 3.5839 | 90 | 1.0x | 1350 |  |
| `ViscoElasticTransformer::_update_outlier_transformation` | 0.2 | 1.5394 | 90 | 1.0x | 1350 |  |
| `ViscoElasticTransformer::_update_smoothing_weights` | 0.0 | 1.8620 | 3 | 1.0x | 45 |  |

*UNRELIABLE = ms/call < 5× timer overhead (25 ns/scope = 0.000025 ms). Excluded from recommendations.*

## Scaling Exponents (log-log fit across 1K / 10K / 100K)

| Label | Mode | k (exponent) | Complexity bucket | Note |
|---|---|---|---|---|
| `RigidRegistration::update` | rigid | 1.22 | k≈1.16 (n·log n) |  |
| `CorrespondenceFilter::update` | rigid | 1.24 | k≈1.16 (n·log n) |  |
| `SymmetricCorrespondenceFilter::update` | rigid | 1.24 | k≈1.16 (n·log n) |  |
| `NeighbourFinder::update` | rigid | 1.27 | k≈1.16 (n·log n) |  |
| `InlierDetector::update` | rigid | 1.00 | k≈1.0 (linear) |  |
| `ViscoElasticTransformer::_update_elastically` | nonrigid | 1.03 | k≈1.0 (linear) |  |
| `ViscoElasticTransformer::_update_outlier_transformation` | nonrigid | 1.06 | k≈1.0 (linear) |  |
| `SymmetricCorrespondenceFilter::update` | nonrigid | 1.23 | k≈1.16 (n·log n) |  |
| `CorrespondenceFilter::update` | nonrigid | 1.23 | k≈1.16 (n·log n) |  |
| `ViscoElasticTransformer::_update_smoothing_weights` | nonrigid | 1.05 | k≈1.0 (linear) |  |
| `NonrigidRegistration::update` | nonrigid | 1.15 | k≈1.0 (linear) |  |
| `NeighbourFinder::update` | nonrigid | 1.26 | k≈1.16 (n·log n) |  |
| `InlierDetector::update` | nonrigid | 1.00 | k≈1.0 (linear) |  |
| `ViscoElasticTransformer::_apply_transformation` | nonrigid | 1.05 | k≈1.0 (linear) |  |
| `ViscoElasticTransformer::_update_viscously` | nonrigid | 1.03 | k≈1.0 (linear) |  |
| `ViscoElasticTransformer::_update_elastically` | pyramid | 1.03 | k≈1.0 (linear) |  |
| `ViscoElasticTransformer::_update_outlier_transformation` | pyramid | 0.96 | k≈1.0 (linear) |  |
| `PyramidNonrigidRegistration::update` | pyramid | 1.15 | k≈1.0 (linear) |  |
| `SymmetricCorrespondenceFilter::update` | pyramid | 1.21 | k≈1.16 (n·log n) |  |
| `ViscoElasticTransformer::_update_smoothing_weights` | pyramid | 1.02 | k≈1.0 (linear) |  |
| `CorrespondenceFilter::update` | pyramid | 1.21 | k≈1.16 (n·log n) |  |
| `NonrigidRegistration::update` | pyramid | 1.16 | k≈1.0 (linear) |  |
| `NeighbourFinder::update` | pyramid | 1.23 | k≈1.16 (n·log n) |  |
| `NonrigidRegistration::update/layer1` | pyramid | 1.17 | k≈1.0 (linear) |  |
| `NonrigidRegistration::update/layer2` | pyramid | 1.17 | k≈1.0 (linear) |  |
| `NonrigidRegistration::update/layer0` | pyramid | 1.13 | k≈1.0 (linear) |  |
| `ViscoElasticTransformer::_apply_transformation` | pyramid | 1.06 | k≈1.0 (linear) |  |
| `InlierDetector::update` | pyramid | 1.00 | k≈1.0 (linear) |  |
| `ViscoElasticTransformer::_update_viscously` | pyramid | 1.03 | k≈1.0 (linear) |  |

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

- **`NonrigidRegistration::update`** (nonrigid, 1K tier): share=41.3%, Amdahl ceiling = 1/(1−0.413) = 1.7×
  - Scaling exponent k: 1.15
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 15 µs (0.0% of label wall time)

- **`RigidRegistration::update`** (rigid, 1K tier): share=35.5%, Amdahl ceiling = 1/(1−0.355) = 1.6×
  - Scaling exponent k: 1.22
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 15 µs (0.0% of label wall time)

- **`SymmetricCorrespondenceFilter::update`** (rigid, 10K tier): share=26.8%, Amdahl ceiling = 1/(1−0.268) = 1.4×
  - Scaling exponent k: 1.24
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 1215 µs (0.1% of label wall time)

- **`SymmetricCorrespondenceFilter::update`** (nonrigid, 100K tier): share=22.1%, Amdahl ceiling = 1/(1−0.221) = 1.3×
  - Scaling exponent k: 1.23
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 3015 µs (0.0% of label wall time)

- **`NonrigidRegistration::update`** (pyramid, 10K tier): share=21.8%, Amdahl ceiling = 1/(1−0.218) = 1.3×
  - Scaling exponent k: 1.16
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 45 µs (0.0% of label wall time)

- **`CorrespondenceFilter::update`** (rigid, 100K tier): share=20.4%, Amdahl ceiling = 1/(1−0.204) = 1.3×
  - Scaling exponent k: 1.24
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 2430 µs (0.0% of label wall time)

- **`NeighbourFinder::update`** (rigid, 100K tier): share=19.2%, Amdahl ceiling = 1/(1−0.192) = 1.2×
  - Scaling exponent k: 1.27
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 2430 µs (0.0% of label wall time)

- **`CorrespondenceFilter::update`** (nonrigid, 100K tier): share=16.0%, Amdahl ceiling = 1/(1−0.160) = 1.2×
  - Scaling exponent k: 1.23
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 6030 µs (0.0% of label wall time)

- **`NeighbourFinder::update`** (nonrigid, 100K tier): share=14.7%, Amdahl ceiling = 1/(1−0.147) = 1.2×
  - Scaling exponent k: 1.26
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 6045 µs (0.0% of label wall time)

- **`SymmetricCorrespondenceFilter::update`** (pyramid, 100K tier): share=13.6%, Amdahl ceiling = 1/(1−0.136) = 1.2×
  - Scaling exponent k: 1.21
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 1365 µs (0.0% of label wall time)

- **`CorrespondenceFilter::update`** (pyramid, 100K tier): share=9.5%, Amdahl ceiling = 1/(1−0.095) = 1.1×
  - Scaling exponent k: 1.21
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 2730 µs (0.0% of label wall time)

- **`NonrigidRegistration::update/layer2`** (pyramid, 10K tier): share=9.3%, Amdahl ceiling = 1/(1−0.093) = 1.1×
  - Scaling exponent k: 1.17
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 15 µs (0.0% of label wall time)

- **`NeighbourFinder::update`** (pyramid, 100K tier): share=8.8%, Amdahl ceiling = 1/(1−0.088) = 1.1×
  - Scaling exponent k: 1.23
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 2820 µs (0.0% of label wall time)

- **`ViscoElasticTransformer::_update_elastically`** (nonrigid, 1K tier): share=8.3%, Amdahl ceiling = 1/(1−0.083) = 1.1×
  - Scaling exponent k: 1.03
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 3000 µs (5.2% of label wall time)

- **`ViscoElasticTransformer::_update_viscously`** (nonrigid, 1K tier): share=8.2%, Amdahl ceiling = 1/(1−0.082) = 1.1×
  - Scaling exponent k: 1.03
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 3000 µs (5.3% of label wall time)

- **`NonrigidRegistration::update/layer1`** (pyramid, 10K tier): share=7.0%, Amdahl ceiling = 1/(1−0.070) = 1.1×
  - Scaling exponent k: 1.17
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 15 µs (0.0% of label wall time)

- **`NonrigidRegistration::update/layer0`** (pyramid, 1K tier): share=6.4%, Amdahl ceiling = 1/(1−0.064) = 1.1×
  - Scaling exponent k: 1.13
  - Suggested tool: OpenMP (ms/call >> fork-join cost)
  - Fork-join cost: 15 µs (0.1% of label wall time)

### Below granularity threshold

Labels with share ≤ 5% or ms/call ≤ 0.03 ms. OpenMP not viable due to
fork-join overhead exceeding potential speedup (reference: beads 9f5, bdt).

`InlierDetector::update`, `ViscoElasticTransformer::_apply_transformation`, `ViscoElasticTransformer::_update_elastically`, `ViscoElasticTransformer::_update_outlier_transformation`, `ViscoElasticTransformer::_update_smoothing_weights`, `ViscoElasticTransformer::_update_viscously`

### Excluded (UNRELIABLE)

Labels with ms/call < 5× timer overhead (25 ns/scope = 0.000025 ms).
Corrected values unreliable — cannot meaningfully distinguish measurement noise from actual cost.

(none)
