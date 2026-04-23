# Fixes: profiling

**Date:** 2026-04-23
**Review passes:** 1 (competing architecture + implementation reviewers, both Opus/Sonnet)
**Change epic:** `meshmonk-modernization-748` (closed)

## Critical

- **`library/src/PyramidNonrigidRegistration.cpp:126-199`** — The `/layerN` `ScopedTimer` currently wraps the entire per-layer body: float-mesh downsample (lines 130–149), target-mesh downsample (lines 152–168), optional `ScaleShifter` (lines 172–179), AND `nonrigidRegistration.update()` (line 193). The bare `NonrigidRegistration::update` label (inner, in `NonrigidRegistration.cpp`) wraps only the ICP loop. So `sum(/layerN) > bare_total` is structurally guaranteed, and the committed report shows `SANITY WARNING` firing in every tier (6.2–9.1% divergence). Both reviewers caught this independently.
  - **Fix:** move the `_tlayer` ScopedTimer construction from the top of the for-loop body down to wrap ONLY the `nonrigidRegistration.update();` call at line 193. After the fix, `sum(/layerN) ≈ bare_total` within normal measurement noise (<1%). The outer `PyramidNonrigidRegistration::update` label then represents pyramid wall-clock; subtracting `sum(/layerN)` from it yields pure pyramid-setup overhead (downsamplers + scale-shift), which is useful signal.

## Important

- **`library/include/meshmonk/profiling.hpp:75-86`** — `Profiler::calibrate(n)` divides by `n` without guarding `n == 0`, triggering UB on `profiling_calibrate(0)` from Python. Add `if (n == 0) return 0;` at the top of the function body.

- **`library/include/meshmonk/profiling.hpp:75-86`** — `calibrate()` records to the real `g_profiler` accumulator then erases `"__calibration__"`. This leaves a window where `profiling_peek()` between `calibrate()` and `reset()` sees a phantom entry. Refactor: construct a local `Profiler tmp;` inside `calibrate()` and run the 1M-iteration loop against `tmp`, not `*this`. Pure function with zero side-effects on the global accumulator.

- **`profiling/run_profile.py:338-346`** — `_complexity_bucket()` maps log-log fit slopes to labels like `"≈O(n^1.0)"` / `"≈O(n^1.16)"` / `"≈O(n^2.0)"`. The fit is on TOTAL wall-clock (not per-call), so these Big-O labels are misleading. Either (a) fit `ms_per_call` vs `n` and report that (more honest), OR (b) rename the column from "Complexity bucket" to "Growth-rate bucket" and drop the `O(n^k)` notation in favor of just `k≈1.0`. Prefer (a) — fits true per-call complexity.

- **`profiling/run_profile.py` (report recommendations section + pyramid share%)** — After the Critical fix above, the outer `PyramidNonrigidRegistration::update` label will still be a strict superset of its `/layerN` children (still includes downsampling + scale-shift work). Treating it as an independent optimization candidate alongside `NonrigidRegistration::update/layerN` double-counts time. Fix: in `_recommendation_list()` (or wherever recommendations are emitted), SKIP `PyramidNonrigidRegistration::update` when emitting pyramid-mode recommendations — the children already cover the relevant optimization targets. Add a report note: "Pyramid outer label is included for reference only; optimization candidates are the inner /layerN and leaf labels."

- **`tests/test_profiling_driver_748_4.py:595-609`** — `test_exits_with_code_2_when_profiling_disabled` re-implements the guard logic inline (`if not mock_mm.profiling_enabled(): sys.exit(2)`) instead of calling `rp.main()`. A bug in `main()`'s guard would not be caught. Fix: patch `parse_args` and `meshmonk` imports, then call `rp.main()`; assert the `SystemExit.code == 2`.

- **`tests/test_profiling_bindings_748_2.py:143-183`** — `test_on_build_peek_does_not_reset` and `test_on_build_dump_resets_accumulator` operate on an EMPTY accumulator (`{} == {}` passes trivially). Peek could secretly reset and the tests would still pass. Fix: populate the accumulator by calling `profiling_calibrate(n=1000)` (note: after the calibrate-side-effect fix above, this approach won't work — instead run a tiny rigid registration with the committed Template/DemoFace meshes). The existing `TestProfilingRoundtrip.test_reset_peek_peek_dump_roundtrip_with_registration` at lines 230–271 already does this correctly but is conditionally skipped without mesh files — change the skip condition to require only a committed 1K mesh (guaranteed present after bead .3).

- **`profiling/run_profile.py:29,713-714`** — `import numpy as np` at module level (line 29) runs BEFORE `os.environ["OMP_NUM_THREADS"] = "1"` inside `main()` (line 713). On systems with OpenBLAS/MKL-backed numpy, the BLAS thread pool initializes with whatever value was in the environment at numpy import time. For the numpy operations used in `generate_report()` (median, polyfit) this is negligible, but it violates the design doc's "controlled single-threaded baseline" intent. Fix: set `os.environ["OMP_NUM_THREADS"] = "1"` and `EIGEN_DONT_PARALLELIZE = "1"` at the TOP of `run_profile.py` (before any imports that might link BLAS), or move `import numpy as np` inside `main()` after the env pinning.

## Minor

- **`profiling/run_profile.py:551,690`** — UNRELIABLE threshold display: `{overhead_ns / 1e6:.4f} ms` renders 34 ns as `0.0000 ms`, useless to a reader. Change format to `:.6f` (yields `0.000034 ms`) or show it in ns: `5 × {overhead_ns} ns = {5 * overhead_ns} ns`.

- **`profiling/run_profile.py:730`** — calibration is a single 1M-iteration measurement. A transient OS scheduling event during that window biases the result. Run `profiling_calibrate(1_000_000)` 3 times, take the median, record in Methodology. Adds ~2s to driver runtime, substantially improves robustness.

- **`docs/perf/hotspot-profile-20260423-eb81379.md`** — report-quality issue: after the Critical and Important fixes above are implemented, regenerate the first hotspot report. The old report stays in the repo as historical reference; the regenerated one gets a new filename (new SHA → new timestamp).

## ADR Updates

No ADR changes needed. All fixes are implementation-level. ADR-004 D1/D2/D3 (FIRM) and D4 (FLEXIBLE) all remain valid as written.

## Discarded

- **Arch-4** (typed `ProfileEntry` nanobind struct): style preference, D4 FLEXIBLE allows either. Nested dict works and existing Python-side code is clean. Not worth the binding-surface churn.
- **Arch-7** (`PUBLIC` vs `PRIVATE` compile_definition scope): consistent with ADR-004 D3 and design doc rationale. Current `PUBLIC` works; `PRIVATE` would also work for current topology. No bug.
- **Arch-8** (UNRELIABLE threshold has no teeth at 34 ns): noted, but the threshold is correct — it just happens no labels are close to it on this host. Not actionable.
- **Arch-9** (no test coverage): the reviewer missed that tests DO exist at `tests/test_profiling_meshes_748_3.py`, `tests/test_profiling_bindings_748_2.py`, `tests/test_profiling_driver_748_4.py`. Confirmed false finding.
- **Arch-DC-1** (four-function surface via `profiling_snapshot(reset=False)`): the peek/dump split is documented in ADR-004 D4 with explicit rationale (footgun prevention). Not worth re-litigating a recent decision.
