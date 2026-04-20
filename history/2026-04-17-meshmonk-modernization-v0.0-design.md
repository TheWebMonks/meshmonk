# Design: MeshMonk v0.0 — Reality Capture & Scaffolding

**Date:** 2026-04-17
**Status:** Draft
**Parent:** [Overview](./2026-04-17-meshmonk-modernization-design.md)
**Decisions:** [ADR-001](../docs/decisions/ADR-001-meshmonk-modernization.md)

---

## Scope

**v0.0 — Reality capture (~1–2 days):** establish trust anchors before refactoring. Work happens on the `meshmonk-modernization` branch in the dedicated worktree `.worktrees/meshmonk-modernization/` (already created from `origin/cli`). All subsequent work — v0.0 through v0.1 — happens on this single branch in that worktree.

**v0.0 contract: no library changes, no code deletions.** This phase is pure scaffolding + asset relocation + harness utilities. Destructive cleanup (MATLAB `.m` file deletion, `library/examples/example.cpp` removal, tutorial bulk cleanup) is deferred to a dedicated v0.1 D3-execution task behind the `pre-modernization` tag.

Concretely v0.0:
1. Cherry-pick the 8-vertex cube fixture from `origin/add-compute-rigid-transform-cli:cli/test_data/rigid_transform/` into `tests/fixtures/cube_rigid/`.
2. Relocate demo assets (`demo/Template.obj`, `demo/demoFace.obj`, `demo/{Demo,Template}FaceLandmarks.csv`, `demo/rigid_transform.txt`) into `data/`. These files are already tracked on the current branch (inherited from `origin/cli`) — this is a `git mv`, not a cross-branch pull.
3. Ingest the already-committed legacy baseline (`demo/rigid_output.obj`, `demo/pyramid_output.obj`, `demo/rigid_transform.txt`) as Tier 3.5 `.npz` goldens under `tests/golden/legacy_baseline/`. Reproducing these from a fresh `origin/cli` build is a reproducibility check (see checklist), not a blocker.
4. Move only tutorial **assets we want to keep** into `data/` (SIMPOP meshes). Do NOT delete the `.m` scripts or the `library/examples/example.cpp` file in v0.0 — those deletions land in v0.1.
5. Write the mesh-compare helper + obj-to-npz converter with pinned signatures and schema.
6. Set up pytest scaffolding (via `pyproject.toml`, not `pytest.ini`, to align with ADR-001 D5).
7. Commit this design doc + ADR, plus pre-commit config (Python-only hooks in v0.0, scoped to new paths), `NOTICE`, `.editorconfig`.

## Ready-to-execute v0.0 kickoff checklist

When execution starts:

- [x] Create modernization worktree: `git worktree add -b meshmonk-modernization .worktrees/meshmonk-modernization origin/cli`. All subsequent v0.0/v0.1 work happens inside that worktree. (Already done 2026-04-17.)
- [ ] **Verify build prerequisites** (already installed on owner's Mac as of 2026-04-17: `/opt/homebrew/bin/cmake`, `/opt/homebrew/bin/ninja`; AppleClang 17 present). Just confirm `cmake --version ≥ 3.26` for v0.1 compatibility. No install step needed.
- [ ] Cherry-pick `cli/test_data/rigid_transform/` from `origin/add-compute-rigid-transform-cli` into `tests/fixtures/cube_rigid/` — **only these six files** (verified present in worktree `.worktrees/cube-fixture/cli/test_data/rigid_transform/`): `input_vertices.txt`, `input_mesh.obj`, `corresponding_features.txt`, `inlier_weights.txt`, `expected_vertices.txt`, `expected_transform.txt`. Exclude `output_mesh_generated.obj` and `output_transform_generated.txt` (runtime artifacts). **Post-copy assertion:** `ls tests/fixtures/cube_rigid/ | wc -l` equals 6 and `ls tests/fixtures/cube_rigid/*_generated.*` finds nothing.
- [ ] **Relocate demo assets** into `data/` via `git mv` (all already tracked on this branch):
  - `demo/Template.obj` → `data/Template.obj`
  - `demo/demoFace.obj` → `data/demoFace.obj`
  - `demo/DemoFaceLandmarks.csv` → `data/DemoFaceLandmarks.csv`
  - `demo/TemplateLandmarks.csv` → `data/TemplateLandmarks.csv`
  - `demo/rigid_transform.txt` → `data/rigid_transform.txt` (Tier 4 soft anchor)
- [ ] **Pre-move hardcoded-path grep** — before `git mv`, run `grep -r -n -E 'demo/(Template|demoFace|DemoFaceLandmarks|TemplateLandmarks|rigid_transform)' --include='*.m' --include='*.cpp' --include='*.hpp' --include='*.py' --include='*.txt' --include='*.md' .` and record the hits in the commit message. Resolution rule: matches inside `matlab/`, `demo/*.m`, and `tutorial/TutorialScripts/*.m` are acceptable breakage (those files are retired under D3 in v0.1); any hit in `library/`, `cli/`, or `README.md` must be patched in the same commit or the move is aborted.
- [ ] **Tutorial asset relocation (v0.0 subset only)**: `git mv tutorial/TutorialData/SimulatedMappedFaces/SIMPOP_*.obj data/simulated_faces/` (drop the 45 `.mtl` files — unused by the algorithm). **Defer** to a v0.1 D3-execution task: deletion of `tutorial/TutorialScripts/*.m`, `demo/*.m`, `matlab/`, `tutorial/TutorialData/Template_Rigidly_Aligned.*`, `tutorial/TutorialData/Targets/demoFace.*`, and the `pre-modernization` tag. **Catch-all for unlisted files**: `tutorial/TutorialData/Landmarks/demoFace.csv`, `tutorial/TutorialData/SimulatedMetadata.xlsx`, `tutorial/TutorialData/Template.{obj,csv}`, `tutorial/TutorialData/Targets/demoFace.{mtl,png}`, `tutorial/TutorialData/Template_Rigidly_Aligned.png`, and `tutorial/TutorialSlides.pptx` stay in place during v0.0 and are triaged in the v0.1 D3 task. Target for v0.0: SIMPOP asset relocation only (≈70 MB relocated, tutorial/ still present).
- [ ] Set up `tests/` root with pytest scaffolding: `tests/conftest.py`, `tests/__init__.py`, and a root `pyproject.toml` with `[tool.pytest.ini_options]` (no separate `pytest.ini`). The `pyproject.toml` at v0.0 carries only pytest config + project metadata; `scikit-build-core` backend is added in v0.1 per ADR-001 D5.
- [ ] Write `tests/utils/mesh_compare.py` — numpy-only, float64 internally, with these pinned signatures:
  ```python
  def rmse(a: np.ndarray, b: np.ndarray) -> float              # same shape (N,3), per-vertex
  def max_vertex_distance(a: np.ndarray, b: np.ndarray) -> float   # L-inf per-vertex, same shape
  def hausdorff_symmetric(a: np.ndarray, b: np.ndarray) -> float   # bidirectional, cKDTree ok
  ```
  Both inputs are cast to `float64` regardless of incoming dtype. **Initial placeholder tolerances** (to be calibrated in v0.1):
  - Rigid transform components: `atol=1e-4`
  - Per-vertex position (mm): `atol=0.01`
  These land as module-level constants `DEFAULT_TRANSFORM_ATOL` and `DEFAULT_VERTEX_ATOL_MM` with a comment noting v0.1 calibration.
- [ ] Write `tests/utils/obj_to_npz.py` — OBJ → npz converter with this **pinned schema** (documented in module docstring):
  ```
  vertices: (N, 3) float64   — required
  faces:    (M, 3) int32     — required
  normals:  (N, 3) float64   — optional, absent key means "not computed"
  ```
- [ ] Sanity-check that `demo/rigid_output.obj` and `demo/pyramid_output.obj` actually differ before treating them as independent goldens — `cmp` should exit non-zero and mean per-vertex displacement between them should be > 0.1 mm. (Cheap guard against a copy-paste error in the commit that produced them.)
- [ ] **Ingest legacy baseline (Tier 3.5, fast path)**: convert the already-committed `demo/rigid_output.obj`, `demo/pyramid_output.obj`, `data/rigid_transform.txt` → `tests/golden/legacy_baseline/*.npz` via `tests/utils/obj_to_npz.py`. No build required. Also relocate the source OBJs from `demo/` into `tests/golden/legacy_baseline/source_obj/` (kept as provenance) via `git mv`.
- [ ] **Reproducibility check (Tier 3.5, optional-but-recommended)** — divergence is information, not a blocker. Prechecks first: `git -C .worktrees/cli rev-parse HEAD` resolves, and `.worktrees/cli/CMakeLists.txt` exists. Then build `origin/cli` out-of-tree to keep the read-only worktree clean — `cmake -S .worktrees/cli -B /tmp/meshmonk-cli-build && cmake --build /tmp/meshmonk-cli-build --parallel` — then run `./meshmonk_cli rigid_reg "$PWD/.worktrees/cli/demo/Template.obj" "$PWD/.worktrees/cli/demo/demoFace.obj" /tmp/rigid.obj --transform_output /tmp/rigid.txt` and the analogous `pyramid_reg` (all outputs written to `/tmp/`, **never into `.worktrees/cli/demo/`** — the paths in `origin/cli:README.md` write into `demo/` and would overwrite the committed baseline). Compare `/tmp/*.obj` against the committed outputs via `mesh_compare.py`. OpenMesh "complex vertex" / "complex edge" warnings on the demo meshes are expected per `origin/cli:README.md` and do not indicate failure. **Record an explicit outcome label** in the commit message: one of `PASSED`, `DIVERGED <rmse_mm>`, `BUILD_FAILED`, `SKIPPED`. **Fallback if build fails on current Xcode**: run MATLAB demos (`test_rigid_registration.m`, `test_pyramid_registration.m`) — same library via MEX — and compare to the committed outputs. If MATLAB is unavailable on the owner's Mac, record `SKIPPED` and move on.
- [ ] Write first Tier 1 analytical test as a failing placeholder (synthetic rigid recovery). This test is **intentionally red at v0.0** — it ships as a skeleton for v0.1 to fill in.
- [ ] Add `.pre-commit-config.yaml` with **Python-only hooks in v0.0, scoped to v0.0 paths only**. Hooks: `ruff`, `ruff-format`, `check-yaml`, `end-of-file-fixer`, `trailing-whitespace`. **Every hook must carry a `files:` regex** restricting it to v0.0-owned paths so `pre-commit run --all-files` cannot rewrite legacy C++/MATLAB/tutorial files and violate the "no library changes" contract. Use this exact regex for all hooks (shared via YAML anchor or repeated inline):
  ```yaml
  files: '^(data/|tests/|docs/|history/|\.pre-commit-config\.yaml|pyproject\.toml|NOTICE|\.editorconfig|README\.md|CLAUDE\.md|AGENTS\.md)'
  ```
  **Defer** `clang-format-16` and `cmake-format` to v0.1 — enabling them now would rewrite every untouched `.cpp`/`.hpp`/`CMakeLists.txt` file on first `pre-commit run --all-files`, violating the "no library changes" contract. **Verification before commit**: `pre-commit run --all-files` must touch **zero** legacy files (`library/`, `matlab/`, `tutorial/`, `cli/`, `CMakeLists.txt`, `demo/`). If it touches any, the `files:` regex is wrong — fix before committing.
- [ ] Add `.editorconfig` (standard: LF, UTF-8, trim trailing whitespace, final newline, 4-space Python, 2-space YAML/JSON). Note: `.editorconfig` is passive (editor-side), it does not rewrite the repo on install — unlike pre-commit hooks.
- [ ] Add `NOTICE` file at repo root with this exact content:
  ```
  MeshMonk
  Copyright 2017 The Web Monks
  Copyright 2026 jsnyde0

  This product includes software developed by The Web Monks
  (https://github.com/TheWebMonks/meshmonk) prior to the 2026
  Python-first modernization, which was carried out by jsnyde0
  (https://github.com/jsnyde0/meshmonk). Both contributions are
  licensed under the Apache License, Version 2.0.
  ```
- [ ] Commit this design doc and ADR-001 (already in tree; ensure staged).
- [ ] **Tag `v0.0`** on the `meshmonk-modernization` branch once the checklist is green. Tag annotation: `v0.0: scaffolding + reality capture. Test suite ships with an intentionally-failing Tier 1 placeholder (pytest exit code 1 expected at this tag); first green Tier 1 lands in v0.1.` No PR to `cli` at this milestone — v0.0 is scaffolding, not shippable; the first PR back to `cli` is at the v0.1 milestone.
- [x] **Repo transfer**: completed 2026-04-17 — `TheWebMonks/meshmonk` → `jsnyde0/meshmonk`. Local `origin` remote updated; old URL continues to redirect. PyPI namespace claim still TODO.
