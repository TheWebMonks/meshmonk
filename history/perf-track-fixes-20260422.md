# Perf-track review fixes — 2026-04-22

Scope: review of what shipped to `main` for epic `meshmonk-modernization-03i`
(benchmark harness infrastructure + ADR-002 decision record). Stashed
optimization code was intentionally not reviewed per the task brief.

## Auto-resolved (folded into artifacts)

1. `benchmarks/README.md` mesh-size section — previous wording described the
   pre-03i.1 broken downsampler semantics (`ratio = target_n / current_n,
   clamped to 0.5 minimum step`). Updated to match actual `_harness.py`
   behaviour (`ratio_to_remove = 1.0 - target_n / current_n`, clamped to
   **max** 0.5 per step).
2. `benchmarks/README.md` quickstart — removed misleading `-k "5K"` filter
   on the baseline-capture command (no tier uses a "5K" label, so the
   filter would select zero scenarios). Left the note about the 5K-vs-7K
   spec mismatch further down the file as the design record.
3. `docs/decisions/ADR-002-perf-track.md` D4 Outcome 9f5 — stash-index
   reference updated. The original text claimed `stash@{0}`; that was true
   when 9f5 closed but shifted to `stash@{1}` when the bdt pre-triage
   instrumentation was subsequently stashed. Corrected and added a
   durable-identifier note (prefer grepping the stash label over trusting
   the index).

## Deferred — require judgment

### Finding 1: 9f5 Outcome table "Speedup @ best threads" column is arithmetically off for non-7K rows

**Why it's an issue**
The column header says "best threads" (faster of OMP=14 / OMP=4) but the
reported numbers for nonrigid@1K/3K appear to reflect **worst** threads:

| Row | OMP=14 speedup | OMP=4 speedup | Reported in ADR |
|---|---|---|---|
| nonrigid@1K | 725/746 = 0.97x | 725/772 = 0.94x | 0.93-0.94x |
| nonrigid@3K | 847/887 = 0.96x | 847/911 = 0.93x | 0.93-0.94x |

The reported range 0.93-0.94x matches the OMP=4 (slower) column only, not
"best threads". The conclusion (below recalibrated 1.3x gate; guard
violated) is unchanged, but the column header or the numbers are
inconsistent. A reader reconstructing the analysis will find the numbers
don't divide cleanly.

**Recommended fix**
In `docs/decisions/ADR-002-perf-track.md` line 172 table, either
(a) relabel the "Speedup @ best threads" column as "Speedup (range across
OMP=4 / OMP=14)" and show both endpoints, or
(b) recompute the column to use the faster post number for each row
consistently.

**Severity**: Minor — factual, but doesn't change the gate verdict.

---

### Finding 2: 9f5 Outcome table missing rigid@1K and rigid@3K rows

**Why it's an issue**
The recalibrated guard clause (line 162) reads "no regression elsewhere
— rigid and pyramid wall-clock at 1K/3K/7K must not increase (+/-2% noise
band)". The Outcome table only lists `rigid@7K`, omitting rigid@1K and
rigid@3K. A future auditor cannot verify whether rigid@1K/3K regressed
without re-running; the evidence for "guard violated on every non-7K
benchmark" is implicit, not shown.

**Recommended fix**
Add rigid@1K and rigid@3K rows to the table at
`docs/decisions/ADR-002-perf-track.md:172-178`, or note explicitly that
they were not re-run (with reason).

**Severity**: Minor — audit-trail completeness.

---

### Finding 3: bdt Outcome "Nonrigid total" column deviates from committed baseline.json

**Why it's an issue**
The bdt table cites nonrigid totals of 757/895/1112 ms at 1K/3K/7K. The
committed `benchmarks/baseline.json` reports 722.7/843.0/1071.2 ms for
the same scenarios. The ~30-50 ms gap is almost certainly the overhead of
the `std::chrono::steady_clock` instrumentation applied during pre-triage
(timing probes inside three hot loops), not a measurement error. But the
table column simply says "Nonrigid total (6-round median, ms)" without
flagging that these numbers are from an *instrumented* build.

A cold reader comparing the bdt table against `baseline.json` will see a
mismatch and waste time chasing it.

**Recommended fix**
At `docs/decisions/ADR-002-perf-track.md:215`, add a footnote to the
"Nonrigid total" column header: "Instrumented-binary wall-clock;
`std::chrono::steady_clock` probes add ~30-50 ms vs the un-instrumented
baseline in `benchmarks/baseline.json`." The **shares** (the column that
matters for the Amdahl conclusion) are unaffected because both numerator
and denominator come from the same instrumented run.

**Severity**: Minor — factual clarity, does not affect the decision.

---

### Finding 4: 9f5 "Deferred for future revisit" assumes the parallelisable share grows at 50K

**Why it's an issue**
Line 196: "The gate to clear at that point is the original >=2x target,
not the 7K-tier recalibration."

The pre-triage Amdahl ceiling at the 7K-tier parallelisable share of 37%
is 1/(1-0.37) ~= 1.59x. For the original 2x gate to be *theoretically*
reachable, the share must rise to >=50% at 50K. That's plausible for a
k-NN loop that scales worse than O(N) with mesh size, but it's not
proven — the deferred clause treats 2x as reachable without caveating
that assumption.

A future engineer unstashing 9f5 against a 50K mesh might measure, say,
42% parallel share and find 2x still unreachable by Amdahl — then
re-close the bead, frustrated that the deferred clause set an
expectation that couldn't be met.

**Recommended fix**
At `docs/decisions/ADR-002-perf-track.md:196`, replace "The gate to clear
at that point is the original >=2x target" with "The gate to clear at
that point is the original >=2x target, **conditional on re-measuring
that `NeighbourFinder::update` consumes >=50% of nonrigid wall-clock at
50K** — otherwise Amdahl still bounds the attainable speedup below 2x
and the gate must be recalibrated again."

**Severity**: Minor — sets clearer expectations for a future reopening,
but does not affect the current decision.

---

### Finding 5: No consolidated "lessons learned" summary

**Why it's an issue**
The three D4 Outcome sections each contain durable structural insights
(fork-join overhead at fine-grained invocation; Amdahl ceiling from
measured parallel share; allocation-hoist wins shrink with mesh size
because compute grows faster than allocation). A reader who opens
ADR-002 in 6 months looking for "what did we learn from the perf track"
has to read and synthesize all three Outcome sections. The cross-bead
pattern — **pre-triage before implementation, honest measurement, close
cleanly when the gate is missed** — is the most valuable legacy of the
track and deserves a top-level summary pointer.

This is a scope judgment; the content is there, just scattered.

**Recommended fix**
Add a short "D4 Track Summary (2026-04-22)" block immediately after the
three D4 Outcome sections but before D5, listing:
- Three beads attempted, zero merged.
- Common structural reasons (cite the Outcomes for detail).
- The pre-triage-first pattern as a template for future perf beads.
- The harness itself is the durable deliverable.

Alternatively, a one-paragraph note at the top of the ADR pointing
readers to the three Outcome sections in reading order would suffice
with less structural change.

**Severity**: Minor — story coherence for future readers.

---

## No fixes required elsewhere

- Harness integrity: `_harness.py::downsample_to` correctly implements
  `ratio_to_remove = 1.0 - target_n / current_n` capped at 0.5 max per
  step. Docstring matches behaviour.
- `conftest.py` xdist guard is hard-assert at import time with a clear
  error message — fails fast as intended.
- OMP_NUM_THREADS pinning is effective; the Eigen::setNbThreads(1) caveat
  is accurately documented (no Python binding exists; env-var form of
  EIGEN_DONT_PARALLELIZE is a no-op at runtime, confirmed against
  `vendor/eigen-3.4.0/Eigen/Core:65`).
- All three bench files use the canonical `load_full_mesh` /
  `downsample_to` helpers from `_harness.py`; no duplicated logic.
- pedantic() config is honest: rounds=5, warmup_rounds=1, consistent
  across rigid/nonrigid/pyramid. `baseline.json` confirms 5 rounds and
  median aggregation.
- `baseline.json` machine context is fully captured (aarch64, 14 cores,
  Apple, 3.11.15, OrbStack Linux). Reproducible.
- 9f5 Outcome baseline numbers match `baseline.json` within ~4 ms (same
  harness, re-run rounding). Internally consistent.
- `InlierDetector::_smooth_inlier_weights` is indeed dead code —
  verified: defined at `library/src/InlierDetector.cpp:84` but never
  called from `InlierDetector::update()` or elsewhere. ADR's dead-code
  claim is correct.
- Cross-linking: README references ADR-002, ADR references README.
- Stash labels are descriptive enough that a stranger could re-measure
  9f5 / bdt pre-triage from `git stash list` output.
