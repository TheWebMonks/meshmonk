# Fixes: neighbourfinder-profiling
Date: 2026-04-23
Review passes: 1

Both reviewers converged on the same three structural issues. Triage below.

## Critical

- **docs/perf/neighbourfinder-breakdown-20260423-807f67c.md:13, 25-30** —
  Committed report has `Repo SHA: unknown` and all six mesh file SHA-256 entries as "file not found". This violates ADR-004 §Methodology's reproducibility contract (must record host, compiler, build flags, git SHA, mesh file SHA-256). A future reader cannot verify which commit or which mesh generated these numbers. Root cause is the `run_local.py` shim bypassing `run_profile.py`'s `/workspace`-bound provenance lookups without substituting its own. Fix flows from Important #1 below: once `run_profile.py` is parameterized on `REPO_ROOT`, rerun and regenerate the report with correct provenance. Old file replaced in place (same commit SHA in filename).

## Important

- **profiling/run_profile.py:140, 448, 451, 463, 824** —
  Five hardcoded `/workspace` paths prevent the documented `python -m profiling.run_profile ...` invocation from working on any host without the `/workspace` symlink. Parameterize on `REPO_ROOT = Path(__file__).parent.parent.resolve()`: derive `DATA = REPO_ROOT / "data"`, `cwd=REPO_ROOT` for the two `git rev-parse` calls in `generate_report`, mesh SHA lookups under `REPO_ROOT / "data"`, and `out_dir = REPO_ROOT / "docs" / "perf"` in `main`. This makes the driver work identically in CI (where `REPO_ROOT == /workspace`) and on dev hosts.

- **profiling/run_local.py (entire file, 191 lines)** —
  Delete. Once `run_profile.py` is parameterized (above), the shim is unnecessary. It currently duplicates argument parsing, the calibration block, the measurement loop, the scaling-fit call, and the output-path logic, while silently diverging from `run_profile.py`'s label-sanity check. Keeping it guarantees drift on any future `run_profile.py` change. Root-cause fix is the parameterization; no reason to keep a second driver.

- **library/src/NeighbourFinder.hpp (3 of 4 ScopedTimer scopes)** —
  Remove `NeighbourFinder::buffer_alloc`, `NeighbourFinder::query_setup`, `NeighbourFinder::result_copy`. Keep only `NeighbourFinder::tree_query`. The report already establishes all three are below the 5× noise floor (UNRELIABLE per ADR-004 D2). They contribute ~3.9s of overhead per 100K run (120.9M timer constructions × 32ns) that contaminates every future `MESHMONK_PROFILING=ON` run at every scale — not just this investigation. The investigation conclusion (`tree_query` dominates) is already captured in the report and committed; keeping instrumentation that cannot measure itself is pure ongoing cost. Retain `tree_query` so the follow-up neighbour-caching bead has a direct before/after signal on the dominant cost bucket. Keep the D2 buffer hoist — it's a positive structural change regardless.

## Minor

- **library/src/NeighbourFinder.hpp (comment on hoisted buffers)** —
  Update the "resize/reinit per-iteration buffers" comment to reflect that `resize()` is a no-op after iteration 0 (vectors stay at target capacity). One-line fix. Avoids a documentation trap for future readers who might try to interpret the hoisted code as per-iteration allocation.

- **docs/perf/neighbourfinder-breakdown-20260423-807f67c.md (regenerate)** —
  After the driver fixes and scope removal above, rerun at 100K nonrigid with the unified driver. The new report will have (a) correct git SHA and mesh SHA-256 provenance, (b) lower instrumentation overhead (40.3M scopes instead of 161.2M), and (c) a cleaner `NeighbourFinder::update` vs `NeighbourFinder::tree_query` share comparison. Filename uses the new commit's short SHA.

## ADR Updates

- **docs/decisions/ADR-006-neighbourfinder-profiling.md (new)** —
  Promote `history/2026-04-23-neighbourfinder-profiling-design.md` to ADR-006 with an Outcome section. Rationale: project convention is ADRs for FIRM investigation decisions (ADR-005 followed this pattern for an investigation that was never implemented). Record:
  - D1 (FIRM, retained): Nested `{}` scope brackets for each cost bucket
  - D2 (FIRM, retained): Hoist buffer declarations above the loop
  - D3 (FIRM, retained): Run at 100K nonrigid only
  - D4 (FIRM, new): Sub-function ScopedTimer labels inside per-vertex loops must be removed after the investigation that added them, OR gated behind a separate verbose flag. Rationale: the `MESHMONK_PROFILING=ON` build is a shared artifact used by all profiling investigations; per-vertex instrumentation adds seconds of overhead that contaminates unrelated measurements.
  - Outcome: `tree_query` dominates; bead for neighbour caching (#1) is the recommended next step per the 5jf decision tree.

  Add Parent link to ADR-004 (profiling infrastructure) and Related link to ADR-005 (superseded OpenMP ADR).

## Discarded

- **Report filename date vs body date mismatch** (Arch #7) — cosmetic; resolves as side-effect of single-driver parameterization (one `datetime.date.today()` call).
- **Dead `--tiers`/`--modes` in run_local.py** (Arch #8) — subsumed by shim deletion.
- **`KNNResultSet` copy-assignment in `buffer_alloc` vs `query_setup` labeling mismatch** (Impl #1) — subsumed by scope removal; both buckets are being deleted.
- **`resize()` no-op semantic mismatch** (Arch #4, Impl #5) — the `buffer_alloc` bucket is being removed, which moots the bucket-intent mismatch. Comment fix (Minor above) handles what remains.
- **Report recommendation confidence** (Impl #4) — reviewers agreed the "Medium confidence" rating is appropriate and not overclaimed. No change needed.

## Pipeline note

Per ADR-063 D4 (FIRM), there is NO re-review after this fix pass. Human catches any remainder.
