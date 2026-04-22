# MeshMonk Benchmark Harness

pytest-benchmark harness for the MeshMonk performance track (ADR-002).

## Quick start

```bash
# Install dev dependencies (includes pytest-benchmark)
uv sync --all-extras

# Smoke run (minimal rounds, verbose)
OMP_NUM_THREADS=1 uv run pytest benchmarks/ --benchmark-only --benchmark-min-rounds=2 -v

# Full baseline capture (all tiers — committed as baseline.json).
# Note: tiers implemented are 1K/3K/7K; there is no 5K label. Use no -k
# filter (or filter on actual tier labels like "1K or 3K or 7K").
OMP_NUM_THREADS=1 uv run pytest benchmarks/ --benchmark-only \
    --benchmark-json=benchmarks/baseline.json \
    --benchmark-min-rounds=5
# Review autonomous fold-in 2026-04-22: previous example used -k "5K"
# which selects zero scenarios since no tier uses that label.

# Compare against baseline (simulate what nightly CI runs)
OMP_NUM_THREADS=1 uv run pytest benchmarks/ --benchmark-only \
    --benchmark-json=new.json \
    --benchmark-compare=benchmarks/baseline.json \
    --benchmark-compare-fail=median:20%
```

## Mesh-size parameterisation

### Decision: Cap at 7K (tier-reduction approach)

`data/Template.obj` has **7,160 vertices**. The original design doc assumed
≥50K vertices, which is not the case. Two approaches were considered:

| Approach | Pros | Cons |
|---|---|---|
| **Cap at 7K (chosen)** | Simple, no extra code, reproducible | Top tier is smaller than ideal |
| Edge-split upsampling to ≥50K | Faithful to design intent | Requires a deterministic subdivision implementation; adds complexity |

The cap-at-7K approach was chosen because:
- Random vertex duplication is forbidden (benchmark inputs must be
  reproducible).
- A correct deterministic edge-split implementation would take longer to
  build than the harness itself.
- The harness goal is regression detection, not absolute throughput numbers
  at industrial scale. 7K is sufficient to characterise the hot paths.

**Tiers:**
- **1K** (~1,000 vertices): obtained by iterative downsampling from the full mesh.
- **3K** (~3,000 vertices): obtained by iterative downsampling.
- **7K** (7,160 vertices): full Template.obj, used as-is.

Downsampling uses `meshmonk.downsample_mesh` with deterministic ratios.
Note: `meshmonk.downsample_mesh(ratio)` interprets `ratio` as the fraction
of vertices to **remove** (i.e. number of edge collapses / current vertex
count), not the fraction to keep. The harness computes
`ratio_to_remove = 1.0 - target_n / current_n` clamped to a **maximum of
0.5 per step** (never remove more than half the mesh in a single pass) and
iterates until the target is reached or decimation stalls. This matches
the implementation in `benchmarks/_harness.py::downsample_to` and the
ratio-semantics fix from bead 03i.1. The ratio is deterministic given the
same input mesh, so benchmark inputs are reproducible across machines.
<!-- Review autonomous fold-in 2026-04-22: previous wording described the
     pre-03i.1 broken semantics ("ratio = target_n / current_n, clamped to
     0.5 minimum step"). Updated to match actual _harness.py behaviour. -->


### The committed baseline (baseline.json)

Only the **5K** scenario label is used in nightly CI (see AC6 in the bead).
Since the mesh tiers are 1K/3K/7K, the "5K" label in `pytest -k "5K"`
selects nothing. The actual baseline committed here captures the **3K** and
**7K** tiers by running without `-k` filtering.

> **Note for future maintainers:** if you want to align the tier labels with
> the bead spec (which says "5K"), either add a 5K tier or rename the 3K
> tier to "5K". For now the baseline reflects the actual tiers implemented.

## Thread control and baseline isolation

### OMP_NUM_THREADS=1

All baseline runs pin `OMP_NUM_THREADS=1` so the baseline captures strictly
single-threaded MeshMonk behaviour. The `conftest.py` `pin_threads` fixture
does this automatically for every session started inside `benchmarks/`.

Post-optimisation benchmarks (after 9f5, bdt land) should be run **without**
`OMP_NUM_THREADS=1` so the reported speedup is the incremental OpenMP
contribution.

### Eigen::setNbThreads(1)

`Eigen::setNbThreads(1)` is the correct runtime call to disable Eigen's own
threading (`EIGEN_DONT_PARALLELIZE` is a compile-time macro; the env-var
form is a no-op at runtime). No Python binding for `setNbThreads` exists in
meshmonk yet. When added, `conftest.py` should be updated to call it. In the
meantime, `OMP_NUM_THREADS=1` also constrains Eigen's OpenMP-backed paths.

## xdist oversubscription guard

Running benchmarks under `pytest-xdist` (`-n auto`) multiplies the OMP
thread budget across workers, producing noisy and non-reproducible results.
`conftest.py` hard-asserts `PYTEST_XDIST_WORKER not in os.environ` and
raises immediately if xdist is active. Never use `-n` when running benchmarks.

## Batch-workload oversubscription footgun

If you call `meshmonk.rigid_register()` (or nonrigid/pyramid) from multiple
Python threads concurrently (e.g. `ThreadPoolExecutor` batching registrations),
each Python thread spawns `OMP_NUM_THREADS` OpenMP threads internally. On an
8-core machine with 8 Python workers and OMP default, this creates 64 threads
for 8 cores — heavy oversubscription.

**Fix:** set `OMP_NUM_THREADS=1` and parallelise at the Python level instead.
MeshMonk is designed for intra-call parallelism, not batch workloads from
multiple concurrent callers.

## Baseline refresh cadence

`benchmarks/baseline.json` must be refreshed:
- On every **minor release** (version bump x.Y.0), OR
- When the **GitHub Actions runner image** changes (`ubuntu-latest` rolls
  ~annually), whichever comes first.

Stale baselines turn the regression gate into noise. The nightly workflow
includes a `refresh_baseline` `workflow_dispatch` input: set it to `true` to
re-capture and open a PR updating `baseline.json`. Do NOT update baseline.json
manually or auto-merge on green — every baseline update requires a PR review
(ADR D6).

## CI regression tolerances

| Registration type | Tolerance | Reason |
|---|---|---|
| rigid | 20% | Deterministic iteration count |
| nonrigid | 20% | Deterministic iteration count |
| pyramid | 30% | Annealing-driven iteration count is noisier on CI runners |

Note: `--benchmark-compare-fail=median:20%` does not support per-scenario
tolerances. The CI workflow uses a post-step script (`check_pyramid_tolerance.py`
in CI) to apply the 30% threshold to pyramid scenarios separately.

## Files

| File | Purpose |
|---|---|
| `conftest.py` | xdist guard, OMP thread pinning, session config |
| `bench_rigid.py` | Rigid registration benchmarks (1K/3K/7K) |
| `bench_nonrigid.py` | Nonrigid registration benchmarks (1K/3K/7K) |
| `bench_pyramid.py` | Pyramid registration benchmarks (1K/3K/7K) |
| `baseline.json` | Committed baseline (single-threaded, OMP_NUM_THREADS=1) |
| `README.md` | This file |
