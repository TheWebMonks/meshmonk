# Fixes: meshmonk-modernization-v0.0
Date: 2026-04-17
Review passes: 2

## Critical

- **`.gitignore`** — Remove the legacy `data/` ignore line (line ~53 from the upstream WebMonks repo, labeled "Exlucde data folder"). It silently drops any NEW file added under `data/` — a footgun for Tier 3.5/4 golden assets. Currently-tracked files survive only because they were `git mv`'d. Fix: delete the `data/` line entirely. Verify with `git check-ignore -v data/__probe.txt` returning nothing.

- **`pyproject.toml`** — Missing `[build-system]` table; orphan `[tool.setuptools]` forces the implicit setuptools backend (visible via the `meshmonk.egg-info/` that had to be `.gitignore`'d). v0.0 is metadata-only, so the cleanest fix is: remove `[tool.setuptools] packages = []` entirely AND add an explicit `[build-system] requires = ["setuptools>=77"]`, `build-backend = "setuptools.build_meta"` so the backend choice is intentional at v0.0. (v0.1 swaps to scikit-build-core per ADR-001 D5.)

## Important

- **`tests/utils/obj_to_npz.py:38-39`** — Meshio stores OBJ vertex normals under key `"obj:vn"`, not `"Normals"`. The current branch is dead code; any OBJ with normals silently loses them. Fix: change `"Normals"` → `"obj:vn"` in both the key check (line 38) and the `point_data[...]` read (line 39).

- **`tests/utils/obj_to_npz.py`** — `convert()` silently drops non-triangle cells (quads, polygons) when triangles are also present. The existing `ValueError("no triangles")` only fires when zero triangles exist. Fix: after extracting triangles, compare total cell count; if any were dropped, `raise ValueError(f"non-triangle cells encountered ({dropped} dropped); retriangulate upstream")`. This matches the pinned "triangles only" contract.

- **`tests/utils/obj_to_npz.py:24,33`** — Arrays are not forced C-contiguous. Downstream v0.1 Eigen::Map consumers will hit hidden copies (contradicts ADR-001 D6 zero-copy intent). Fix: use `np.ascontiguousarray(mesh.points, dtype=np.float64)` for vertices, `np.ascontiguousarray(np.concatenate(...), dtype=np.int32)` for faces, same for normals. Document C-contiguity in the module docstring schema.

- **`tests/utils/mesh_compare.py`** — Three different failure modes on empty input: `rmse` → NaN (silent), `max_vertex_distance` → `ValueError` from numpy, `hausdorff_symmetric` → `ValueError` from cKDTree. Fix: add a shared guard at each function entry: `if a.shape[0] == 0 or b.shape[0] == 0: raise ValueError("empty input array")`. Also validate `a.ndim == 2 and a.shape[1] == 3` so a `(N,2)` or `(N,4)` array fails loud instead of returning a wrong number.

- **`tests/utils/test_obj_to_npz.py` — `test_roundtrip_values`** — Asserts one vertex is present but never checks face indices. The MINIMAL_OBJ fixture has deterministic faces `[[0,1,2],[0,1,3],[0,2,3],[1,2,3]]` (0-based after conversion). Fix: add `assert np.array_equal(sorted(data["faces"].tolist()), [[0,1,2],[0,1,3],[0,2,3],[1,2,3]])` (or equivalent set-membership). Also add an end-to-end normals roundtrip test to catch regressions of the key-name bug above.

- **`tests/utils/test_mesh_compare.py`** — `max_vertex_distance` and `hausdorff_symmetric` only have identity/symmetry tests — a stub returning a constant would pass. Fix: add `test_max_vertex_distance_known` using e.g. `a=[[0,0,0],[1,0,0]]`, `b=[[0,0,0],[2,0,0]]` → expected 1.0; add `test_hausdorff_known_value` using `a=[[0,0,0],[1,0,0]]`, `b=[[0.5,0,0]]` → expected 0.5.

- **`pyproject.toml`** — Dev extras (numpy, scipy, meshio, pytest) unpinned. `uv.lock` pins them but raw pip users get resolver roulette, which defeats the v0.0 reality-capture premise. Fix: add minimum floors matching the lock: `numpy>=1.26`, `scipy>=1.11`, `meshio>=5.3`, `pytest>=8.0`.

- **Tag `v0.0` collides with legacy `v0.0.1`–`v0.0.6`** — Seven pre-existing tags from the conda/MATLAB era share the same prefix; `v0.0.6` sorts above `v0.0` semver-wise. Fix: add a section to `docs/v0.0-release-notes.md` and ADR-001 explicitly calling out the naming collision ("v0.0–v0.3 refers to the modernization track; legacy v0.0.1–v0.0.6 are archived MATLAB/conda releases and should not be re-tagged"). Alternative (heavier): move v0.0 → `modernization-v0.0`. Prefer doc-only fix now since the tag is already pushed.

## Minor

- **`tests/test_rigid_analytical.py`** — Use `@pytest.mark.xfail(raises=NotImplementedError, strict=True)` instead of an unmarked failing test. Keeps CI green (when it exists) while still enforcing the stub contract. Update `docs/v0.0-release-notes.md` to note "1 xfailed" instead of "1 failed".

- **`tests/fixtures/cube_rigid/README.md`** — Missing. `legacy_baseline/` has one; this dir should match the convention. Fix: add a README with provenance (upstream branch/SHA: `cli/test_data/rigid_transform/` on `origin/add-compute-rigid-transform-cli`), the 6-file manifest, and the purpose (Tier 2 analytical fixture).

- **`tests/golden/legacy_baseline/README.md`** — Add a `## Reproducibility check (v0.0)` section recording: outcome `BUILD_FAILED` (acceptable per design), root cause (`cxxopts` rejects `f` suffix on float literals in `cli/cli.cpp:139,141,154,156`), pointer to `docs/v0.0-reproducibility-check.md` for full details, pointer to tracking bead `meshmonk-modernization-eh1`.

- **`pyproject.toml`** — Missing `description` and `readme` fields under `[project]`. Fix: add `description = "Python-first 3D mesh registration library"`, `readme = "README.md"`.

- **`tests/conftest.py`** — Empty except for a one-line comment. Fix: delete the file entirely; pytest doesn't need it, and empty-with-a-promise files invite drift. (v0.1 can add it back when it actually has fixtures.)

## ADR Updates

- **ADR-001** — Add a short note under the v0.0 section: "Tag naming collides with legacy `v0.0.1`–`v0.0.6` (conda/MATLAB era). The modernization track reuses the prefix intentionally; legacy tags are archived and will not be re-published." No decision change; clarification only.

## Discarded

- Pre-commit allowlist comment explaining v0.1 extension path (design doc already sets the expectation; comment would rot).
- `tests/stubs/rigid.py` relocation / docstring handoff line (trivial; v0.1 swap is a mechanical import rename).
- Legacy `Makefile` deprecation banner (v0.0 contract forbids touching it; v0.1 removes or wraps it).
- `docs/developing.md` (premature; v0.1 will have an actual install story with scikit-build-core).
- `.editorconfig` `tab_width = 4` cosmetic (no impact with `indent_style = space`).
- `issues.jsonl` post-commit churn (pre-existing beads-hook hygiene issue, not introduced by v0.0; file a separate repo-level bead if it keeps recurring).
- No-CI finding (Critical-sounding but ADR-001 already lists it; CI belongs with v0.1 when there's a green Tier 1 to gate).
