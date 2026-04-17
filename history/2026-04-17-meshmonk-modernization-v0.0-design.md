# Design: MeshMonk v0.0 — Reality Capture & Scaffolding

**Date:** 2026-04-17
**Status:** Draft
**Parent:** [Overview](./2026-04-17-meshmonk-modernization-design.md)
**Decisions:** [ADR-001](../docs/decisions/ADR-001-meshmonk-modernization.md)

---

## Scope

**v0.0 — Reality capture (~1–2 days):** establish trust anchors before refactoring. Create the `meshmonk-modernization` branch from `origin/cli` and check it out in the dedicated worktree `.worktrees/meshmonk-modernization/` (see Branch strategy). All subsequent work — v0.0 through v0.1 — happens on this single branch in that worktree. Cherry-pick the 8-vertex cube fixture from `origin/add-compute-rigid-transform-cli:cli/test_data/rigid_transform/` to `tests/fixtures/cube_rigid/` before that branch is retired. Move current demo assets (`demo/Template.obj`, `demo/demoFace.obj`, `demo/{Demo,Template}FaceLandmarks.csv`) to `data/`, plus `origin/cli:demo/rigid_transform.txt` (verified present on the `cli` branch, not in current checkout). **Capture legacy baseline:** the `cli` branch already has committed outputs — `origin/cli:demo/rigid_output.obj`, `origin/cli:demo/pyramid_output.obj`, `origin/cli:demo/rigid_transform.txt` — from a prior run of `meshmonk_cli` on `Template.obj ↔ demoFace.obj`. These become the v0.0 Tier 3.5 reference, converted to `.npz` via `tests/utils/obj_to_npz.py` and stored under `tests/golden/legacy_baseline/`. Reproducing them from a fresh `origin/cli` build on the owner's Mac is a **reproducibility check** (see checklist), not a blocker for v0.0 — `cmake` is not currently installed locally. Write the mesh-compare helper. Set up pytest scaffolding. Commit this design doc + ADR. No library changes.

## Ready-to-execute v0.0 kickoff checklist

When execution starts:

- [ ] Create modernization worktree: `git worktree add -b meshmonk-modernization .worktrees/meshmonk-modernization origin/cli`. All subsequent v0.0/v0.1 work happens inside that worktree.
- [ ] Install build prerequisites on owner's Mac: `brew install cmake ninja` (verified absent as of 2026-04-17). AppleClang 17 already present and satisfies C++14 requirement for `origin/cli`.
- [ ] Cherry-pick `cli/test_data/rigid_transform/` from `origin/add-compute-rigid-transform-cli` into `tests/fixtures/cube_rigid/` — **only these files** (verified present in worktree `.worktrees/cube-fixture/cli/test_data/rigid_transform/`): `input_vertices.txt`, `input_mesh.obj`, `corresponding_features.txt`, `inlier_weights.txt`, `expected_vertices.txt`, `expected_transform.txt`. Exclude `output_mesh_generated.obj` and `output_transform_generated.txt` (runtime artifacts from the old CLI, not fixture inputs).
- [ ] Move `demo/Template.obj`, `demo/demoFace.obj`, `demo/{Demo,Template}FaceLandmarks.csv` → `data/`. Also pull `origin/cli:demo/rigid_transform.txt` into `data/` (present on `cli`, absent in current checkout — verified 2026-04-17) as the Tier 4 soft anchor.
- [ ] **Clean up `tutorial/` (99 MB)**: (a) retire `tutorial/TutorialScripts/*.m` alongside `demo/*.m` and `matlab/` under D3 (tag `pre-modernization` before deletion); (b) move `tutorial/TutorialData/SimulatedMappedFaces/SIMPOP_*.obj` → `data/simulated_faces/` (drop the 45 `.mtl` files — unused by the algorithm); (c) delete `tutorial/TutorialData/Template_Rigidly_Aligned.*` and `tutorial/TutorialData/Targets/demoFace.*` (duplicates of demo/ assets); (d) defer `TutorialSlides.pptx` until a licensed-clean alternative or owner-held slides exist. Target: ≥80 MB removed from working tree.
- [ ] Delete root `example.cpp` and any branch-local example scratchpads (Python CLI is the canonical reference; C++ consumers link against the installed library).
- [ ] Set up `tests/` root with pytest scaffolding (`conftest.py`, empty `__init__.py`, minimal `pyproject.toml` or `pytest.ini`)
- [ ] Write `tests/utils/mesh_compare.py` — Hausdorff / RMSE / per-vertex distance helpers (~20 lines, numpy only)
- [ ] Write `tests/utils/obj_to_npz.py` — OBJ → npz converter for Tier 3.5 legacy-baseline ingestion
- [ ] Sanity-check that `origin/cli:demo/rigid_output.obj` and `demo/pyramid_output.obj` actually differ before treating them as independent goldens — `cmp` should exit non-zero and mean per-vertex displacement between them should be > 0.1 mm. (Cheap guard against a copy-paste error in the commit that produced them.)
- [ ] **Ingest legacy baseline (Tier 3.5, fast path)**: convert already-committed `origin/cli:demo/rigid_output.obj`, `demo/pyramid_output.obj`, `demo/rigid_transform.txt` → `tests/golden/legacy_baseline/*.npz` via `tests/utils/obj_to_npz.py`. No build required.
- [ ] **Reproducibility check (Tier 3.5, optional-but-recommended)**: build `origin/cli` out-of-tree to keep the read-only worktree clean — `cmake -S .worktrees/cli -B /tmp/meshmonk-cli-build && cmake --build /tmp/meshmonk-cli-build --parallel` — then run `./meshmonk_cli rigid_reg "$PWD/.worktrees/cli/demo/Template.obj" "$PWD/.worktrees/cli/demo/demoFace.obj" /tmp/rigid.obj --transform_output /tmp/rigid.txt` and the analogous `pyramid_reg` (all outputs written to `/tmp/`, **never into `.worktrees/cli/demo/`** — the paths in `origin/cli:README.md` write into `demo/` and would overwrite the committed baseline). Compare `/tmp/*.obj` against the committed outputs via `mesh_compare.py`. OpenMesh "complex vertex" / "complex edge" warnings on the demo meshes are expected per `origin/cli:README.md` and do not indicate failure. Divergence is information — document, do not block. **Fallback if build fails on current Xcode**: run MATLAB demos (`test_rigid_registration.m`, `test_pyramid_registration.m`) — same library via MEX — and compare to the committed outputs.
- [ ] Write first Tier 1 analytical test as a failing placeholder (synthetic rigid recovery)
- [ ] Add `.pre-commit-config.yaml` (ruff, clang-format-16, cmake-format), `.editorconfig`, `NOTICE` file (WebMonks 2017 + jsnyde0 2026 dual-copyright)
- [ ] Commit this design doc and ADR-001
- [ ] Tag `v0.0` on the `meshmonk-modernization` branch once the checklist above is green (no PR to `cli` at this milestone — v0.0 is scaffolding, not shippable; the first PR back to `cli` is at the v0.1 milestone)
- [x] **Repo transfer**: completed 2026-04-17 — `TheWebMonks/meshmonk` → `jsnyde0/meshmonk`. Local `origin` remote updated; old URL continues to redirect. PyPI namespace claim still TODO.
