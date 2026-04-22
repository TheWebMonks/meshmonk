# ADR-002: Performance Track — Harness-First, OpenMP-Parallel

**Status:** Accepted
**Date:** 2026-04-21
**Design:** [Perf track design doc](../../history/2026-04-21-perf-track-design.md)
**Parent:** [ADR-001 meshmonk-modernization](./ADR-001-meshmonk-modernization.md)
**Related:** Beads `meshmonk-modernization-44m`, `-9f5`, `-6a5`, `-bdt`

## Context

The registration hot paths (`NeighbourFinder::update`, `ViscoElasticTransformer::_update_*`, `InlierDetector::_smooth_inlier_weights`, `NonrigidRegistration`) have never been benchmarked or parallelized. Four related beads were filed 2026-04-21 after a triage scan. The work needs a coherent sequencing + invariants framework, not four independent patches.

OpenMP is not currently linked — `CMakeLists.txt:44-45` explicitly notes this. Adding it is a real but minor toolchain change, and it's the natural choice given Eigen already uses OpenMP when available and the targeted loops are embarrassingly parallel.

This ADR records the decisions that shape the perf epic. Algorithmic identity and the public API surface are not up for discussion here — they're covered by ADR-001.

## Decisions

### D1: Measure before optimize — benchmark harness is a strict prerequisite

**Firmness: FIRM**

No performance optimization merges to main without (a) a baseline recorded on the `benchmarks/` harness, and (b) a measured win over that baseline demonstrated in the PR description.

**Rationale:**

The codebase has no benchmark history. Perf claims without measurement accumulate as untested folklore. The harness (`meshmonk-modernization-44m`) must land first to give every subsequent change a real number to beat. This also sets up a regression gate that catches accidental slowdowns in unrelated changes — cheap insurance.

The harness uses `pytest-benchmark` (added to `[project.optional-dependencies].dev` in `pyproject.toml` by `44m` — it is NOT currently in the dev stack despite being low-ceremony to adopt) with scenarios covering rigid/nonrigid/pyramid × (1K, 5K, 50K) vertices. A curated subset runs nightly on CI with a 20% regression tolerance for rigid/nonrigid and a 30% tolerance for pyramid (variable iteration count from annealing makes it noisier). Baselines run with `OMP_NUM_THREADS=1` and `EIGEN_DONT_PARALLELIZE=1` pinned so the reported speedup is strictly the incremental effect of MeshMonk-level pragmas, not contaminated by Eigen's opportunistic threading.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **pytest-benchmark harness first (chosen)** | • Low ceremony<br>• Measurable wins<br>• Regression gate | • ~1 day upfront before any perf work |
| Optimize first, benchmark later | • Faster to first patch | • Perf claims unverifiable<br>• No regression protection<br>• Invites drift |
| asv (airspeed velocity) | • Richer history UI | • Heavier setup<br>• Overkill for 4-scenario matrix |
| google-benchmark for C++ only | • Native C++ timing | • Doesn't exercise the Python boundary (where users care) |

**What would invalidate this:**

- A perf bug materializes that's visibly broken (>10× slowdown) and needs a same-day hotfix — then we ship the fix with a post-hoc benchmark rather than blocking on harness-first.

---

### D2: OpenMP for parallelization — not std::execution, not TBB

**Firmness: FIRM**

Adopt OpenMP (`#pragma omp parallel for`) for the parallelizable loops in `NeighbourFinder::update`, `ViscoElasticTransformer::_update_*`, and `InlierDetector::_smooth_inlier_weights`. Add `find_package(OpenMP REQUIRED)` + `target_link_libraries(meshmonk_lib PUBLIC OpenMP::OpenMP_CXX)` to `library/CMakeLists.txt`.

**Rationale:**

- Eigen already links OpenMP when available (see `vendor/eigen-3.4.0/Eigen/src/Core/products/Parallelizer.h`). Adding it at the meshmonk level doesn't introduce a novel dependency concept — it formalizes what Eigen was already probing for.
- OpenMP is shipped with every target compiler (GCC, Clang via libomp, MSVC) and requires zero runtime install on Linux/Windows. macOS needs `libomp` via Homebrew; cibuildwheel already handles this.
- The targeted loops are embarrassingly parallel — no reductions, no shared-mutable-state, output writes at distinct row indices. OpenMP's pragma-based model is the right level of abstraction; `std::execution` adds C++17 parallel algorithms complexity without benefit.
- TBB would add a runtime dependency to the wheel. Not worth it for three loops.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **OpenMP (chosen)** | • Pragma-based, minimal code churn<br>• Already probed by Eigen<br>• No runtime dep on Linux/Windows | • macOS needs libomp via brew<br>• Non-deterministic parallel reductions |
| `std::execution::par` (C++17 parallel algorithms) | • Standard C++ | • Loops don't map cleanly to parallel algorithms (they're not `for_each` over a range with a simple predicate)<br>• GCC requires TBB backend anyway |
| Intel TBB | • Mature, work-stealing | • Runtime dep ships in wheel<br>• Higher binding complexity |
| Hand-rolled `std::thread` pool | • No external dep | • Reinvent scheduling, reinvent bugs |
| Keep sequential | • No toolchain change | • Leaves 4-8× on the table |

**What would invalidate this:**

- Apple Silicon `libomp` proves unreliable in cibuildwheel (verify during 9f5).
- A future GPU backend is introduced — then CUDA/HIP parallelism replaces OpenMP for those paths, but OpenMP stays for CPU fallback.

---

### D3: Tolerance-based numeric invariance check, not bitwise

**Firmness: FIRM**

Every perf change must pass the existing e2e test suite (`tests/test_registration_e2e.py`) with its current tolerances. Bitwise-equal output is **not** required — parallel float reductions are non-associative and bit-reproducibility across thread counts is neither achievable nor worth pursuing.

**Rationale:**

`6a5` (cache temporaries, no parallelism) should produce bitwise-identical output since the operations are unchanged — and D4 enforces this specifically (see D4 criterion 5). `9f5` and `bdt` introduce OpenMP and may produce output that differs by a few ULPs on some entries.

The existing e2e suite in `tests/test_registration_e2e.py` uses scalar aggregate thresholds (`step.mean() < 10.0`, `move.mean() > 10.0`, `diff.mean() < 5.0`, `step.max() < 50.0`) — these catch gross breakage but are too coarse to detect subtle loop-variable races introduced by incorrect parallelization. The **parallel-vs-sequential parity test** described in the design doc (a dedicated test run of the affected registration path under `OMP_NUM_THREADS=1` vs default threads, asserting output equivalence within ~1e-5 on aligned vertex coordinates) is the real numeric guard for `9f5` and `bdt`, and is a **required artifact** of those beads — not a "nice to have".

If a user needs deterministic output, `OMP_NUM_THREADS=1` restores sequential behavior. This is documented in `benchmarks/README.md` and a future `docs/gotchas.md`.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **Tolerance-based (chosen)** | • Matches physical reality of parallel fp<br>• No artificial constraints on implementation | • Users wanting bitwise determinism must set `OMP_NUM_THREADS=1` |
| Bitwise-identical required | • Fully deterministic | • Infeasible with parallel reductions without perf cost<br>• Forces sequential or deterministic-reduction schemes |
| Deterministic reduction (e.g. tree-reduce with fixed associativity) | • Deterministic & parallel | • Adds complexity to every loop<br>• Slower than naive parallel sum |

**What would invalidate this:**

- A downstream user requires bitwise reproducibility for regulatory reasons — would need a separate deterministic mode, probably via compile-time flag.
- E2e tolerances are tightened below the parallel-fp noise floor.

---

### D4: Per-optimization merge gate — measurable win required

**Firmness: FIRM**

Each of `9f5`, `6a5`, `bdt` merges only when the PR description includes:
1. Baseline wall-clock from `benchmarks/baseline.json` for affected scenarios.
2. Post-change wall-clock on the same harness.
3. E2e tests green (with existing aggregate-threshold tolerances).
4. Speedup ≥ a bead-specific threshold: `9f5` ≥ 2× on 5K-nonrigid, `6a5` ≥ 5% on 5K-pyramid (raised from an initial 3% based on prior GitHub Actions run-to-run variance data — 3% falls inside CI noise), `bdt` ≥ 1.5× on 5K-nonrigid. The `6a5` revised expected range (3–8% after the F1/F2 inner-loop-only scope boundary) means this bead may fall below its threshold; that is acceptable per this ADR's own escape hatch (below) — the PR closes and we re-triage with real numbers. Refresh the floor after 44m provides a noise characterization.
5. **For `6a5` only:** bitwise-identical e2e output vs baseline. `6a5` changes memory layout, not operations — output must not shift. A diff of registration vertex arrays vs a committed reference is cheap insurance against accidental reordering during the hoist. `9f5` and `bdt` are exempt from bitwise equality (parallel float reductions are non-associative) and instead must pass the parallel-vs-sequential parity test (see D3).

**Rationale:**

Without per-bead thresholds, we risk merging a "parallel" change that's actually slower due to thread-launch overhead on small inputs. The thresholds are set conservatively to give each change room to underperform expectations while still being net-positive.

If a bead can't meet its threshold, the PR is closed (not merged with caveats) and the bead re-triaged — maybe the hot path isn't actually hot, maybe the parallelization is wrong, maybe the harness scenario is miscalibrated. That's information we want surfaced cleanly, not buried in merged code.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **Per-bead threshold gate (chosen)** | • Forces honest measurement<br>• Closes ineffective changes cleanly | • Slight friction on borderline wins |
| Merge any measured win | • Lower friction | • ≤5% "wins" often lost in CI noise; accumulating them wastes review budget |
| No threshold, just require green tests | • Simplest | • Invites "parallel but not actually faster" merges |

**What would invalidate this:**

- The harness proves too noisy to discriminate at the 3% level — then `6a5`'s threshold raises to 5-10%.

### D4 Outcome (2026-04-22): 6a5 closed below gate

**Measurement summary:** 1.3–4.0% speedup across 9 scenarios (rigid/nonrigid/pyramid × 1K/3K/7K tiers) on Template.obj (7160 vertices max).

**Key result:** Pyramid@7K = 3.0%. The 7K tier is the closest achievable proxy for the ADR's 5K-pyramid gate target due to Template.obj tier-snapping; the 50K scenario from D1 requires a mesh we do not have.

**Gate outcome: MISS** — 3.0% vs ≥5% threshold. No scenario reached the gate.

**Scaling analysis:** Allocation-hoist wins shrink with mesh size. Per-iteration allocation cost is O(N); compute is O(N log N) or worse. As N grows, allocation's share of wall-clock falls. A 50K synthetic upsample would not be expected to rescue the gate — the scaling direction is unfavorable.

**Note on prior measurement:** A prior harness run measured 5.75% pyramid@7K, but that was against a broken harness — the downsampler had a ratio-semantics defect fixed in bead 03i.1. That number is not comparable to the post-03i.1 results above.

**Action per D4 protocol:** PR closed (code not merged), bead closed. Code preserved in `git stash` labeled "6a5 hoist code — closed below D4 gate 2026-04-22 (see ADR-002 D4 Outcome)" rather than destroyed, in case future benchmark coverage (e.g. a real 50K mesh, or different workload profile) reopens the question.

**Retained artifacts:** `tests/test_memory_layout_regression.py` + `tests/golden/memory_layout_reference.npy` kept as reusable insurance for future memory-layout changes to InlierDetector, NonrigidRegistration, and ViscoElasticTransformer. Originally captured as the 6a5 parity test; reframed as a generic regression guard (no 6a5-specific framing).

### D4 Gate Recalibration (2026-04-22): 9f5 threshold lowered to ≥1.3× on nonrigid@7K

**Original gate:** `9f5 ≥ 2× on 5K-nonrigid` (criterion 4 above).

**Recalibrated gate:** `9f5 ≥ 1.3× on nonrigid@7K, with no regression elsewhere`.

**Measurement summary (pre-triage, 2026-04-22, post-03i.1 harness):** The parallelizable k-NN path (`NeighbourFinder::update`) accounts for 26–37% of nonrigid wall-clock at the 1K/3K/7K tiers (37% at 7K). By Amdahl's law, parallelising that fraction alone caps total speedup at `1 / (1 − 0.37) ≈ 1.59×` even with a hypothetically free parallel KD-tree search. A realistic parallel KD-tree gets 3–4× on that fraction, giving an expected total speedup of ~1.36–1.45× at 7K — well below the original 2× threshold but comfortably above the recalibrated 1.3× floor.

**Why 2× is unreachable for this workload mix:** The original 2× target implicitly assumed the parallelizable share was ≥50%. Post-harness measurement showed it is not. The remaining non-parallel share (Eigen linear solves, host-side transform application, per-iteration setup) is the dominant cost, and is out of 9f5's scope.

**Mesh-size substitution:** 7K rather than 5K, for the same Template.obj tier-snapping reason documented in the 6a5 Outcome above (we do not have a 5K-exact mesh, and we lack a real 50K mesh for the D1 high-end tier).

**Additional guard clause:** "no regression elsewhere" — rigid and pyramid wall-clock at 1K/3K/7K must not increase (±2% noise band). Parallel thread-launch overhead on smaller inputs (1K) is the failure mode we want to catch.

**If the recalibrated gate is also missed:** Close per the standard D4 protocol (stash + document + close bead). The recalibration is not a gate softening to let the work through — it is an honest correction of an analytical assumption that turned out not to match the code.

---

### D4 Outcome (2026-04-22): 9f5 closed below recalibrated gate

**Measured result against the recalibrated `≥1.3× on nonrigid@7K, no regression elsewhere` gate:**

| Benchmark | Baseline (OMP=1) | Post (OMP=14) | Post (OMP=4) | Speedup @ best threads |
|---|---|---|---|---|
| **nonrigid@7K (gate)** | 1069 ms | 1111 ms | 1105 ms | **0.97×** (MISS vs ≥1.3×) |
| nonrigid@1K / 3K | 725 / 847 ms | 746 / 887 ms | 772 / 911 ms | 0.93–0.94× (regression, guard violated) |
| pyramid@1K / 3K / 7K | 419 / 457 / 545 ms | 433 / 485 / 565 ms | 430 / 468 / 553 ms | 0.97–0.99× (regression, guard violated) |
| rigid@7K | 1117 ms | — | 1162 ms | 0.96× (regression, guard violated) |

Both gate criteria failed: gate threshold unmet (0.97× < 1.3×), and the no-regression guard clause triggered on every non-7K-nonrigid benchmark.

**Implementation was correct:** parallel-vs-sequential parity test (D3 required artifact) produced bitwise-identical aligned vertex coordinates under `OMP_NUM_THREADS=1` vs `OMP_NUM_THREADS=4` (max-abs-diff = 0.00e+00). All existing e2e and memory-layout tests passed. The OpenMP pragma on `NeighbourFinder::update` compiled, linked, and executed — it was simply slower than the serial baseline.

**Structural reason the optimization did not pay off at benchmark-reachable mesh sizes:**

1. `NeighbourFinder::update` is invoked ~42 times per nonrigid run (one per iteration, for both correspondence directions).
2. Per invocation, the parallelizable workload is ~7000 query points at the largest benchmark tier, split across 14 cores = ~500 pts/thread. Even at 4 threads, it is ~1750 pts/thread.
3. OpenMP parallel-region entry/exit (fork-join + implicit barrier + per-thread stack alloc for the moved-inside temporaries) has fixed overhead per invocation. At this granularity that overhead exceeds the compute savings — especially at the 1K tier, where per-thread workload is <100 points.
4. Reducing thread count from 14 → 4 recovered ~6 ms on a 1069 ms baseline — confirming oversubscription is not the root cause; fork-join cost at the invocation rate is.

The pre-triage Amdahl ceiling of 1.59× assumed idealized parallelization efficiency (zero fork-join cost, linear scaling inside the region). Reality on this workload mix: the k-NN path consumes 26–37% of wall-clock, but parallelising it introduced more overhead than it saved on the benchmark meshes.

**Action per D4 protocol:** Bead `meshmonk-modernization-9f5` closed (code not merged). Working-tree changes preserved in `git stash` labeled "9f5 OpenMP NeighbourFinder — closed below D4 recalibrated gate 2026-04-22 (see ADR-002 D4 Outcome 9f5)". Stash index at time of writing was `stash@{0}`; it subsequently shifted to `stash@{1}` when the bdt pre-triage instrumentation was stashed (see bdt Outcome below for the final index map). The stash **label** is the durable identifier — prefer `git stash list` + grep on the label over a fixed index. Code kept rather than destroyed in case future benchmark coverage (a real 50K mesh, or a workload with larger N per NeighbourFinder invocation) reopens the question.
<!-- Review autonomous fold-in 2026-04-22: corrected stash-index reference
     to match post-bdt stash list ordering; pointed readers to the stash
     label for a durable lookup. -->


**Implications for `bdt`:** The same fork-join overhead concern applies to any OpenMP pragma on an inner loop of a multi-iteration registration. `bdt` targets the nonrigid smoothing/viscoelastic-transform loop. Before implementation, `bdt` must pre-triage the smoothing-loop wall-clock share against the fixed harness (same pattern as 9f5's 2026-04-22 pre-triage) and verify the per-invocation parallelisable work is large enough to amortize OpenMP fork-join cost at current mesh tiers. If pre-triage shows the same pattern as 9f5 (share <50%, fine-grained invocation), `bdt` should be closed or deferred without implementation.

**Deferred for future revisit:** When the D1 high-end tier gains a real 50K-vertex benchmark mesh, re-run 9f5 (unstash, rebuild, measure). At 50K points split across 14 threads ≈ 3600 pts/thread, compute-per-invocation should clear the fork-join break-even. The gate to clear at that point is the original `≥2×` target, not the 7K-tier recalibration.

---

### D4 Outcome (2026-04-22): bdt closed without implementation (pre-triage failed)

Per the 9f5 Outcome's `bdt` clause, a pre-triage measurement was run before any implementation. Result: **all four target loops combined account for 1.7–7.7% of nonrigid wall-clock** — well below the 50% threshold. Closing `meshmonk-modernization-bdt` without implementation.

**Pre-triage methodology:** `std::chrono::steady_clock` instrumentation around the inner `for (size_t i = 0; i < _numElements; ...)` loop in each of the four target functions:

1. `ViscoElasticTransformer::_update_viscously` (library/src/ViscoElasticTransformer.cpp:167)
2. `ViscoElasticTransformer::_update_elastically` (library/src/ViscoElasticTransformer.cpp:214)
3. `ViscoElasticTransformer::_update_outlier_transformation` (library/src/ViscoElasticTransformer.cpp:258)
4. `InlierDetector::_smooth_inlier_weights` (library/src/InlierDetector.cpp:95)

Instrumented binary run against `benchmarks/bench_nonrigid.py` at 1K/3K/7K tiers with `OMP_NUM_THREADS=1`.

**Wall-clock share per tier (inner loop only, cumulative across all ~20 ICP iterations × annealed viscous/elastic/outlier subiterations):**

| Tier | _numElements | Nonrigid total (6-round median, ms) | Viscous (share) | Elastic (share) | Outlier (share) | Smooth inlier (share) | **Combined 4-loop share** |
|------|--------------|-------------------------------------|-----------------|-----------------|-----------------|----------------------|--------------------------|
| 1K   | 1 000        | 757 ms/call                         | 6.3 ms (0.83%)  | 6.0 ms (0.79%)  | 0.28 ms (0.04%) | 0.0 (dead code)       | **1.66%**                |
| 3K   | 3 000        | 895 ms/call                         | 18.3 ms (2.04%) | 18.4 ms (2.05%) | 0.58 ms (0.06%) | 0.0 (dead code)       | **4.16%**                |
| 7K   | 7 000        | 1112 ms/call                        | 42.4 ms (3.81%) | 42.3 ms (3.80%) | 1.32 ms (0.12%) | 0.0 (dead code)       | **7.73%**                |

**Amdahl ceiling computation (7K, the most favourable tier):**

```
upper_bound_speedup = 1 / (1 − parallelizable_share)
                    = 1 / (1 − 0.0773)
                    = 1.084×
```

Even with infinite cores and zero fork-join cost, `bdt` cannot exceed **1.08× on nonrigid@7K** — nowhere near the `≥1.3×` gate that 9f5 was recalibrated to, let alone the original `≥1.5×` gate bdt was scoped to. **The gate is arithmetically unreachable.**

**Per-invocation work is also fine-grained (same regime that defeated 9f5):**

- Viscous inner loop: 7.7 µs/invocation at 1K, 22 µs at 3K, 52 µs at 7K.
- Elastic inner loop: essentially identical to viscous.
- Outlier inner loop: 0.9–4.4 µs across tiers.

9f5 established ~10–30 µs fork-join overhead on this host. Viscous at 1K runs *below* the overhead floor; 3K is approximately at parity; only 7K is modestly above. Total inner-loop fork-joins per `nonrigid_register` call at default annealing: ≈1800 (≈800 viscous + ≈800 elastic + ≈200 outlier). At ~15 µs overhead × 1800 joins ≈ 27 ms of pure scheduling overhead — larger than the entire 1K target workload (12.6 ms combined viscous+elastic) and comparable to the 3K workload (36.7 ms).

**Structural finding: `InlierDetector::_smooth_inlier_weights` is dead code.** The function is declared in the header and defined in the .cpp, but is *never invoked* from `InlierDetector::update()` or anywhere else. It contributes 0.0% to wall-clock because it never runs. This should be noted but does not change the bdt decision; the other three loops at 7.73% combined already fail the gate.

**Narrower-scope analysis:** Is there any subset of the four loops that could pay off?

| Scope | Max share | Upper-bound speedup | Gate |
|-------|-----------|---------------------|------|
| All 4 loops                         | 7.73% | 1.084× | FAIL |
| Viscous only                        | 3.81% | 1.040× | FAIL |
| Elastic only                        | 3.80% | 1.040× | FAIL |
| Viscous + elastic                   | 7.61% | 1.082× | FAIL |
| Outlier only                        | 0.12% | 1.001× | FAIL |

No narrower scope rescues the gate. Every subset fails by wide margin.

**Action per D4 protocol:** Bead `meshmonk-modernization-bdt` closed (code never written — pre-triage preempted implementation). No stash entry for bdt OpenMP code; only a throwaway stash for the pre-triage `std::chrono` instrumentation remains at `stash@{0}` (9f5 OpenMP shifted to `stash@{1}`, 6a5 hoist to `stash@{2}`).

**Follow-ups flagged (not part of this ADR, separate beads if pursued):**

1. **Dead-code removal for `InlierDetector::_smooth_inlier_weights`** — no caller anywhere in the codebase.
2. **~92% of nonrigid wall-clock at 7K is outside the bdt target set.** That time is spent in: `CorrespondenceFilter::update` (k-NN via NeighbourFinder — already measured in 9f5 pre-triage at 26–37% share), `InlierDetector::update` distance loop, Eigen per-row operations in `_update_viscously`'s force-field setup (lines 141–148, which precede the timed inner loop), and `_apply_transformation`. Future perf work should pre-triage these before assuming OpenMP is the right tool; 9f5 demonstrated that fine-grained parallelism on short-lived invocations does not amortize.

**Deferred for future revisit:** If the D1 high-end tier gains a real 50K-vertex benchmark mesh, re-measure the smoothing shares. At 50K vertices, viscous/elastic per-invocation work would be ~370 µs (7.3× 7K cost) — well above fork-join overhead. But the *share* of total nonrigid time depends on how the correspondence filter scales; the shares could stay low even at 50K. Revisit only after the 50K mesh exists and the overall nonrigid wall-clock profile at that size is known.

---

### D5: `OMP_NUM_THREADS` is the only public thread knob for v0.x

**Firmness: FLEXIBLE**

Thread count is controlled via the standard `OMP_NUM_THREADS` environment variable. No Python-level knob (`n_threads=`, `set_num_threads()`) is added in this track.

**Rationale:**

Users who care about thread control know `OMP_NUM_THREADS`. Adding a Python surface means adding binding code, documentation, conflict resolution with the env var, and testing — none of which is justified until a user actually asks. YAGNI.

If demand materializes (e.g. web-worker contexts where env vars are awkward), we can add `meshmonk.set_num_threads(n)` as a thin wrapper over `omp_set_num_threads()` in a follow-up bead. Backward compatible.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **OMP_NUM_THREADS only (chosen)** | • Zero added surface<br>• Standard convention | • Users unfamiliar with OpenMP may not know to set it |
| Add `meshmonk.set_num_threads()` | • Discoverable | • Unjustified complexity until demand |
| Per-call `n_threads=` kwarg | • Most flexible | • Threads through every register function; API bloat |

**What would invalidate this:**

- A user explicitly requests per-call or library-level thread control (then we add a thin wrapper).
- We end up in a context where env vars don't work (rare for a library; flag if encountered).

---

### D6: Benchmark baseline JSON is committed, PR-gated refresh

**Firmness: FLEXIBLE**

`benchmarks/baseline.json` lives in the repo, captures the current harness results, and only updates via an explicit PR that re-runs the harness and commits the new JSON. No auto-refresh on green builds.

**Rationale:**

A drifting baseline is a useless regression gate. Forcing explicit refresh gives the person updating it a chance to explain *why* the baseline is moving (new hardware, intentional optimization, legitimate regression being accepted).

The baseline is tied to the CI runner profile (a specific GitHub Actions image). Documented in `benchmarks/README.md`. Local runs won't match the baseline — that's fine; local is for iteration, CI is for regression detection.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **Committed JSON, PR-gated refresh (chosen)** | • Baseline is auditable<br>• Forces intentional updates | • Mild PR friction on legitimate refreshes |
| Auto-update on green main | • Zero friction | • Baseline drifts; gate becomes meaningless |
| External benchmark service (asv, Airspeed Velocity) | • History UI, more insight | • Heavier ops burden<br>• Overkill for 4 scenarios |

**What would invalidate this:**

- We outgrow 4-scenario coverage and want per-commit history — then migrate to asv.

---

### D7: No scope creep — algorithmic or API changes belong elsewhere

**Firmness: FIRM**

This track does not change the public API, does not alter algorithms, does not introduce new dependencies beyond OpenMP, and does not address docstring/type-hint/test-harness concerns (covered by separate beads `-hnf`, `-xvu`, `-h31`, `-0bs`, `-ald`, `-d3q`, `-gw1`, `-4x4`).

**Rationale:**

The four perf beads are cohesive precisely because they share the measurement methodology and OpenMP toolchain shift. Bundling unrelated polish would dilute the review surface and blur the merge gate criteria. Each of the other 8 beads can be sent through their own pipeline once this track lands and baseline numbers exist.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **Perf-only scope (chosen)** | • Clean merge gates<br>• Focused review | • Temptation to sneak in small fixes must be resisted |
| Bundle with API polish (hnf/xvu/h31) | • Fewer PRs | • Mixes "measurable perf" with "judgment-call API" review lenses |
| Bundle with C++ test harness (0bs) | • Some overlap (CMake changes) | • 0bs has different acceptance criteria (test coverage, not speedup) |

**What would invalidate this:**

- A discovered bug in the hot paths requires an algorithmic fix before parallelization is safe — then that fix lands first, separately, and this track rebases on it.

## Related

- [Perf track design doc](../../history/2026-04-21-perf-track-design.md) — implementation details
- [ADR-001 meshmonk-modernization](./ADR-001-meshmonk-modernization.md) — parent strategic ADR; D1 (C++20), D5 (CMake + scikit-build-core) both relevant here
- Beads: `meshmonk-modernization-44m`, `-9f5`, `-6a5`, `-bdt`
