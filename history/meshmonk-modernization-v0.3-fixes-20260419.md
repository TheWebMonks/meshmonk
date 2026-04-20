# Fixes: meshmonk-modernization-v0.3
Date: 2026-04-19
Review passes: 2 (architecture + implementation reviewers per pass)

## Critical

- **pyproject.toml** — Missing `[tool.scikit-build] sdist.exclude` config. scikit-build-core includes all `git ls-files` output by default. Without `sdist.exclude = ["tests/", "data/", "history/"]`, the `build-sdist` job's self-check (`grep -E '(tests|data)/'`) will fail, blocking all publish jobs. Fix: add `[tool.scikit-build]` section with `sdist.exclude`.

- **pyproject.toml:32** — `pyyaml` missing from cibuildwheel `test-requires`. `test_integration_v03.py` has `import yaml` at module level. cibuildwheel only installs `["pytest", "scipy", "trimesh"]`, so pytest collection fails with `ModuleNotFoundError`. Fix: add `"pyyaml"` to `test-requires` list. Also add to `test-sdist` job's pip install in `release.yml`.

- **.github/workflows/release.yml:59,93,114-215** — `actions/upload-artifact` pinned with comment `# v7.0.1` and `actions/download-artifact` with `# v8.0.1`. As of 2026, highest published major version for both is v4. The SHAs are likely fabricated and will fail at runtime with "unable to resolve action." Fix: resolve actual SHAs for `actions/upload-artifact@v4` and `actions/download-artifact@v4`, update all references.

## Important

- **meshmonk/cli.py:226** — `urllib.request.urlretrieve` has no timeout and no error handling. Network failures expose raw `URLError` tracebacks. Fix: wrap in `try/except (urllib.error.URLError, OSError)`, emit user-friendly message via `typer.echo`, raise `typer.Exit(code=1)`. Use `urlopen(url, timeout=30)` with manual file write instead of `urlretrieve`.

- **meshmonk/cli.py:221-240** — Partial file leak on download failure. If `urlretrieve` fails mid-download, truncated file remains in `_CACHE_DIR`. Subsequent `_find_demo_meshes()` finds it (only checks `.exists()`) and passes to trimesh, causing confusing parse error. Fix: wrap download in try/except, `dest.unlink(missing_ok=True)` on any exception before re-raising.

- **docs/api/types.md:48-56** — `RigidTransform` docs list `.rotation` (ndarray (3,3)), `.translation` (ndarray (3,)), `.scale` (float) attributes. None exist in the binding — only `.matrix` (4x4 ndarray), `.compose()`, `.apply()`, `.inverse()`. Users following docs get `AttributeError`. Fix: replace with actual API.

- **.github/workflows/release.yml:14-16** — `id-token: write` and `contents: write` granted at workflow level, giving all jobs (including `build-wheels` which runs third-party code) OIDC token access. Fix: set `permissions: {}` at workflow level, move permissions into only `publish-testpypi`, `publish-pypi`, and `create-release` jobs.

- **.github/workflows/release.yml:198-203** — `create-release` depends on `publish-pypi` via `needs`. When `publish-pypi` is skipped (RC/dev tags), `create-release` is also skipped. RC tags never get GitHub Releases. Fix: change `needs` to `[build-wheels, build-sdist]` so create-release runs independently of publish jobs. Add explicit RC/dev filtering to `create-release`'s `if` condition if RC releases should NOT get GitHub Releases.

- **.github/workflows/release.yml:172-176** — `publish-pypi` job condition requires `github.event_name == 'push'`, making the `workflow_dispatch` option `pypi` inert. Fix: either remove `pypi` from workflow_dispatch options (only `testpypi` and `none`), or add workflow_dispatch condition to publish-pypi.

## Minor

- **docs/api/errors.md:22-29** — Example shows only `print(exc)` but not `.code` property — the primary mechanism for programmatic error dispatch (added in v0.2). Fix: add example showing `if exc.code == RegistrationError.DegenerateInput:` pattern.

- **docs/api/types.md:12-19** — Primary `aligned_features` field (N,6) undocumented for result types. Only derived `aligned_vertices` shown. Users who need normals alongside positions won't know `aligned_features` exists. Fix: document `aligned_features` as the primary field.

## ADR Updates

- No ADR changes needed. All findings are implementation issues, not decision revisions.

## Discarded

- **ci.yml not SHA-pinned**: design doc specifies SHA pinning for release.yml only. Not v0.3 scope.
- **Validation boilerplate dup (~69 lines)**: explicitly deferred to v0.4+ per design doc prerequisites table.
- **aarch64 test-skip**: performance concern. Acceptable for first release, can optimize later.
- **GH Pages deployment missing from release.yml**: design doc explicitly says docs site "can ship after the first release."
- **/releases/latest/ vs sha256 conflict**: design choice — sha256 is intentionally None until release, at which point both are set.
- **DEMO_ASSETS sha256 None**: acknowledged prerequisite, can't compute until release exists.
- **floating_flags Pattern A restriction**: pre-existing API design, not v0.3 change.
- **test_demo_runs_rigid_mode monkeypatch**: pre-existing test quality issue.
- **Substring match in rc/dev tag detection**: edge case requiring unusual tag names, low risk.
- **_apply_pyramid_kwargs explicit_kwargs redundant**: pre-existing code, not v0.3 change.
- **from __future__ import annotations**: no functional issue.
- **Rigid pre-alignment normals**: mathematically correct behavior (normals rotate with rigid transform).
- **test_ci_test_command_excludes_silent_logger false-passing**: minor test quality, not blocking.
- **mkdocs renders full ADR**: only 1 ADR exists, acceptable for v0.3.
- **create-release condition matches rc/dev tags**: correctly relies on needs chain to publish-pypi.
