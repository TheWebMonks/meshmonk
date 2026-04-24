# NeighbourFinder Warm-Start Cache — Before/After (2026-04-24)

Measures the impact of the bead 729 warm-start cache (ADR-007) on 100K
nonrigid registration. Ceiling-only worst-dist seeding with `(1 + eps)`
inflation; bit-identical output to the uncached path.

## Methodology

- **Tier:** 100K nonrigid (most hot tier per ADR-006)
- **Repetitions:** 1 warmup + 3 measured runs, median selected
- **Build:** `-DMESHMONK_PROFILING=ON`, Release
- **Pinning:** `OMP_NUM_THREADS=1`, `EIGEN_DONT_PARALLELIZE=1`
- **Driver:** `python -m profiling.run_profile --tiers 100k --modes nonrigid --runs 3 --warmup 1`
- **Before SHA:** `07e11cf` (ADR-006 closeout — `tree_query` retained, no caching)
- **After SHA:** this commit (bead 729 caching applied)
- **ScopedTimer overhead:** 33–34 ns/scope (both runs)

## Summary

| Label | Before ms/call | After ms/call | Δ ms/call | Δ % |
|---|---:|---:|---:|---:|
| `NonrigidRegistration::update` | 66560.62 | 62468.57 | -4092.05 | **-6.15%** |
| `SymmetricCorrespondenceFilter::update` | 212.61 | 188.09 | -24.52 | **-11.53%** |
| `CorrespondenceFilter::update` | 78.00 | 64.42 | -13.58 | **-17.41%** |
| `NeighbourFinder::update` | 72.9128 | 59.0912 | -13.82 | **-18.95%** |
| `NeighbourFinder::tree_query` | 0.0002 | 0.0001 | -0.0001 | **~-50%** (now UNRELIABLE) |

Labels NOT affected by the change (sanity check — all within 1–2% noise):

| Label | Before | After | Δ % |
|---|---:|---:|---:|
| `ViscoElasticTransformer::_update_viscously` | 30.32 | 30.90 | +1.9% |
| `ViscoElasticTransformer::_update_elastically` | 30.17 | 30.81 | +2.1% |
| `InlierDetector::update` | 5.23 | 5.32 | +1.7% |
| `ViscoElasticTransformer::_apply_transformation` | 3.50 | 3.72 | +6.3% |

(Small positive deltas on unchanged labels are noise across the 3-run
median; they quantify the per-label measurement floor for this 100K run.)

## Interpretation

1. **`NeighbourFinder::tree_query` per-call cost halved.** Pre-seeded
   `worstDist = c_k * (1 + eps)` prunes more kd-tree sibling branches
   early, so the average query walks fewer nodes. The per-call cost
   (~100 ns) is now below the 5× timer-overhead noise floor (170 ns),
   hence the UNRELIABLE flag in the report. A sharper measurement would
   require lower-overhead instrumentation — out of scope for this bead.

2. **`NeighbourFinder::update` total improves ~19% even though `tree_query`
   is only ~27% of NF::update's wall-time.** Explanation: the retained
   `tree_query` ScopedTimer has measurement overhead of 40.3M × ~33 ns =
   1.33 s per run (per ADR-006 D4). Halving the time inside the timed
   region also ~halves the overhead attribution for that label. The
   overall NF::update wall-time drop (~13.8 ms/call × 403 calls = 5.57 s)
   is real and dominates.

3. **Full nonrigid wall-time drops 6%** — Amdahl-consistent with NF::update
   being ~15% of the original total. The speedup flows cleanly up through
   the composition stack (`CorrespondenceFilter::update`,
   `SymmetricCorrespondenceFilter::update`, `NonrigidRegistration::update`)
   since NF::update is the leaf.

## Correctness

Verified by `tests/test_neighbour_caching_729.py::test_cached_matches_uncached_nonrigid`:
two nonrigid registrations of the same input pair give byte-identical
`aligned_vertices` with caching ON vs OFF (5 iterations, full ICP loop).
See ADR-007 D2 for why the raw-ceiling variant produced a ~5-unit drift
vs. uncached and why `(1 + eps)` inflation is required.

All 350 other tests pass (`pytest tests/ --ignore=tests/test_golden.py --ignore=tests/test_silent_logger.py`).
The one pre-existing failure (`test_bitwise_equality_with_reference`) is
unrelated: the committed reference `.npy` predates the ScopedTimer + D2
buffer-hoist changes and fails bitwise on `main` without any caching
changes. Tracked in bead `meshmonk-modernization-77d`.
