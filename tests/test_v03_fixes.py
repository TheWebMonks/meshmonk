"""Tests for v0.3 review fixes (bead workspace-eeg.5).

TDD: these tests are written BEFORE the implementation changes.
They cover: C1, C2, I1, I2, I3, I4, I5, M1, M2.
"""
from __future__ import annotations

import re
from pathlib import Path

import yaml

try:
    import tomllib
except ImportError:
    import tomli as tomllib  # type: ignore[no-redef]

ROOT = Path(__file__).parent.parent


# ---------------------------------------------------------------------------
# C1: sdist.exclude in pyproject.toml
# ---------------------------------------------------------------------------


def test_pyproject_has_scikit_build_section():
    """pyproject.toml must have [tool.scikit-build] section."""
    content = tomllib.loads((ROOT / "pyproject.toml").read_text())
    assert "scikit-build" in content.get("tool", {}), (
        "pyproject.toml missing [tool.scikit-build] section"
    )


def test_pyproject_scikit_build_sdist_exclude():
    """[tool.scikit-build] must have sdist.exclude with tests/, data/, history/."""
    content = tomllib.loads((ROOT / "pyproject.toml").read_text())
    skbuild = content["tool"]["scikit-build"]
    sdist = skbuild.get("sdist", {})
    exclude = sdist.get("exclude", [])
    assert "tests/" in exclude, "sdist.exclude must include 'tests/'"
    assert "data/" in exclude, "sdist.exclude must include 'data/'"
    assert "history/" in exclude, "sdist.exclude must include 'history/'"


# ---------------------------------------------------------------------------
# C2: pyyaml in cibuildwheel test-requires
# ---------------------------------------------------------------------------


def test_cibuildwheel_test_requires_has_pyyaml():
    """cibuildwheel test-requires must include pyyaml."""
    content = tomllib.loads((ROOT / "pyproject.toml").read_text())
    test_requires = content["tool"]["cibuildwheel"]["test-requires"]
    assert any("pyyaml" in r.lower() for r in test_requires), (
        f"pyyaml missing from cibuildwheel test-requires: {test_requires}"
    )


def test_release_yml_test_sdist_installs_pyyaml():
    """release.yml test-sdist job pip install step must include pyyaml."""
    content = (ROOT / ".github/workflows/release.yml").read_text()
    # Find the test-sdist section and check pyyaml appears there
    parsed = yaml.safe_load(content)
    test_sdist_job = parsed["jobs"]["test-sdist"]
    steps = test_sdist_job["steps"]
    # Look for a pip install step that includes pyyaml
    pip_steps = [
        s for s in steps
        if s.get("name", "").lower().startswith("install test")
        or ("run" in s and "pip install" in str(s.get("run", "")))
    ]
    found_pyyaml = any("pyyaml" in str(s.get("run", "")).lower() for s in pip_steps)
    assert found_pyyaml, (
        "test-sdist job must have a pip install step that includes pyyaml"
    )


# ---------------------------------------------------------------------------
# C3: upload-artifact/download-artifact must use v4-compatible SHAs
# ---------------------------------------------------------------------------


def test_release_yml_no_fabricated_v7_v8_artifact_actions():
    """release.yml must not reference actions/upload-artifact or download-artifact with v7/v8 comments."""
    content = (ROOT / ".github/workflows/release.yml").read_text()
    assert "v7.0.1" not in content, (
        "release.yml still references fabricated upload-artifact v7.0.1 SHA"
    )
    assert "v8.0.1" not in content, (
        "release.yml still references fabricated download-artifact v8.0.1 SHA"
    )


def test_release_yml_artifact_actions_are_v4():
    """All upload-artifact and download-artifact references must be v4."""
    content = (ROOT / ".github/workflows/release.yml").read_text()
    upload_refs = re.findall(r"actions/upload-artifact@(\S+)", content)
    download_refs = re.findall(r"actions/download-artifact@(\S+)", content)
    for ref in upload_refs:
        # Either a valid 40-char SHA OR v4 floating tag
        sha_valid = re.fullmatch(r"[0-9a-f]{40}", ref) is not None
        tag_valid = ref.startswith("v4")
        assert sha_valid or tag_valid, (
            f"upload-artifact ref is not v4 SHA or tag: @{ref}"
        )
    for ref in download_refs:
        sha_valid = re.fullmatch(r"[0-9a-f]{40}", ref) is not None
        tag_valid = ref.startswith("v4")
        assert sha_valid or tag_valid, (
            f"download-artifact ref is not v4 SHA or tag: @{ref}"
        )


# ---------------------------------------------------------------------------
# I1/I2 (cli.py): timeout + error handling in _demo_download
# ---------------------------------------------------------------------------


def test_cli_download_uses_urlopen_not_urlretrieve():
    """cli.py _demo_download must use urlopen (with timeout) not urlretrieve."""
    content = (ROOT / "meshmonk/cli.py").read_text()
    assert "urlretrieve" not in content, (
        "_demo_download still uses urllib.request.urlretrieve (no timeout)"
    )
    assert "urlopen" in content, (
        "_demo_download must use urllib.request.urlopen with timeout=30"
    )


def test_cli_download_has_timeout():
    """cli.py _demo_download must pass timeout=30 to urlopen."""
    content = (ROOT / "meshmonk/cli.py").read_text()
    assert "timeout=30" in content, (
        "_demo_download missing timeout=30 in urlopen call"
    )


def test_cli_download_catches_urlerror():
    """cli.py _demo_download must catch URLError and OSError."""
    content = (ROOT / "meshmonk/cli.py").read_text()
    assert "URLError" in content, (
        "_demo_download missing URLError in except clause"
    )


def test_cli_download_cleanup_partial_file():
    """cli.py _demo_download must call dest.unlink(missing_ok=True) on exception."""
    content = (ROOT / "meshmonk/cli.py").read_text()
    assert "unlink(missing_ok=True)" in content, (
        "_demo_download must clean up partial file on exception"
    )


def test_cli_download_raises_typer_exit():
    """cli.py _demo_download must raise typer.Exit(code=1) on network error."""
    content = (ROOT / "meshmonk/cli.py").read_text()
    # There should be a raise typer.Exit(code=1) in the error handling path
    # (already exists for SHA mismatch, but must also exist for download errors)
    assert content.count("raise typer.Exit(code=1)") >= 2, (
        "_demo_download must raise typer.Exit(code=1) for both SHA mismatch and network errors"
    )


# ---------------------------------------------------------------------------
# I2: Fix RigidTransform docs
# ---------------------------------------------------------------------------


def test_types_md_rigid_transform_has_matrix_attr():
    """docs/api/types.md RigidTransform must document .matrix attribute."""
    content = (ROOT / "docs/api/types.md").read_text()
    assert "`.matrix`" in content or "| `matrix`" in content, (
        "types.md RigidTransform missing .matrix attribute"
    )


def test_types_md_rigid_transform_no_fake_attrs():
    """docs/api/types.md RigidTransform must not document .rotation, .translation, .scale."""
    content = (ROOT / "docs/api/types.md").read_text()
    # Find the RigidTransform section
    rt_start = content.find("### RigidTransform")
    assert rt_start != -1, "RigidTransform section not found"
    # Find the next section after RigidTransform
    next_section = content.find("\n---", rt_start + 1)
    if next_section == -1:
        next_section = len(content)
    rt_section = content[rt_start:next_section]
    assert "`.rotation`" not in rt_section and "| `rotation`" not in rt_section, (
        "types.md RigidTransform still lists fake .rotation attribute"
    )
    assert "`.translation`" not in rt_section and "| `translation`" not in rt_section, (
        "types.md RigidTransform still lists fake .translation attribute"
    )
    assert "`.scale`" not in rt_section and "| `scale`" not in rt_section, (
        "types.md RigidTransform still lists fake .scale attribute"
    )


def test_types_md_rigid_transform_has_compose_apply_inverse():
    """docs/api/types.md RigidTransform must document .compose(), .apply(), .inverse()."""
    content = (ROOT / "docs/api/types.md").read_text()
    rt_start = content.find("### RigidTransform")
    next_section = content.find("\n---", rt_start + 1)
    if next_section == -1:
        next_section = len(content)
    rt_section = content[rt_start:next_section]
    assert "compose" in rt_section, "RigidTransform missing .compose() method"
    assert "apply" in rt_section, "RigidTransform missing .apply() method"
    assert "inverse" in rt_section, "RigidTransform missing .inverse() method"


# ---------------------------------------------------------------------------
# I3: Per-job permissions in release.yml
# ---------------------------------------------------------------------------


def test_release_yml_workflow_level_permissions_empty():
    """release.yml must have permissions: {} at workflow level."""
    content = (ROOT / ".github/workflows/release.yml").read_text()
    parsed = yaml.safe_load(content)
    workflow_perms = parsed.get("permissions")
    assert workflow_perms == {} or workflow_perms is None, (
        f"Workflow-level permissions must be empty dict or absent, got: {workflow_perms}"
    )
    # Specifically, id-token and contents must NOT be at workflow level
    if workflow_perms:
        assert "id-token" not in workflow_perms, (
            "id-token: write must NOT be at workflow level"
        )
        assert "contents" not in workflow_perms, (
            "contents: write must NOT be at workflow level"
        )


def test_release_yml_publish_testpypi_has_id_token():
    """publish-testpypi job must have id-token: write permission."""
    content = (ROOT / ".github/workflows/release.yml").read_text()
    parsed = yaml.safe_load(content)
    job = parsed["jobs"]["publish-testpypi"]
    perms = job.get("permissions", {})
    assert perms.get("id-token") == "write", (
        "publish-testpypi job missing id-token: write permission"
    )


def test_release_yml_publish_pypi_has_id_token():
    """publish-pypi job must have id-token: write permission."""
    content = (ROOT / ".github/workflows/release.yml").read_text()
    parsed = yaml.safe_load(content)
    job = parsed["jobs"]["publish-pypi"]
    perms = job.get("permissions", {})
    assert perms.get("id-token") == "write", (
        "publish-pypi job missing id-token: write permission"
    )


def test_release_yml_create_release_has_contents_write():
    """create-release job must have contents: write permission."""
    content = (ROOT / ".github/workflows/release.yml").read_text()
    parsed = yaml.safe_load(content)
    job = parsed["jobs"]["create-release"]
    perms = job.get("permissions", {})
    assert perms.get("contents") == "write", (
        "create-release job missing contents: write permission"
    )


def test_release_yml_build_wheels_has_no_sensitive_permissions():
    """build-wheels job must not have id-token: write or contents: write."""
    content = (ROOT / ".github/workflows/release.yml").read_text()
    parsed = yaml.safe_load(content)
    job = parsed["jobs"]["build-wheels"]
    perms = job.get("permissions", {})
    assert perms.get("id-token") != "write", (
        "build-wheels job must NOT have id-token: write"
    )
    assert perms.get("contents") != "write", (
        "build-wheels job must NOT have contents: write"
    )


# ---------------------------------------------------------------------------
# I4: create-release depends on build-wheels + build-sdist (not publish-pypi)
# ---------------------------------------------------------------------------


def test_release_yml_create_release_needs_build_jobs():
    """create-release must depend on build-wheels and build-sdist, not publish-pypi."""
    content = (ROOT / ".github/workflows/release.yml").read_text()
    parsed = yaml.safe_load(content)
    needs = parsed["jobs"]["create-release"]["needs"]
    if isinstance(needs, str):
        needs = [needs]
    assert "build-wheels" in needs, (
        "create-release must need build-wheels"
    )
    assert "build-sdist" in needs, (
        "create-release must need build-sdist"
    )
    assert "publish-pypi" not in needs, (
        "create-release must NOT need publish-pypi (blocks RC releases)"
    )


# ---------------------------------------------------------------------------
# I5: Remove 'pypi' from workflow_dispatch options
# ---------------------------------------------------------------------------


def test_release_yml_no_pypi_dispatch_option():
    """workflow_dispatch publish_target options must not include 'pypi'."""
    content = (ROOT / ".github/workflows/release.yml").read_text()
    parsed = yaml.safe_load(content)
    triggers = parsed.get("on", parsed.get(True, {}))
    wd = triggers.get("workflow_dispatch", {})
    inputs = wd.get("inputs", {})
    options = inputs.get("publish_target", {}).get("options", [])
    assert "pypi" not in options, (
        f"workflow_dispatch options must not include 'pypi' (it's inert): {options}"
    )
    assert "testpypi" in options, "testpypi must remain in workflow_dispatch options"
    assert "none" in options, "none must remain in workflow_dispatch options"


# ---------------------------------------------------------------------------
# M1: .code example in errors.md
# ---------------------------------------------------------------------------


def test_errors_md_has_code_dispatch_example():
    """docs/api/errors.md must show programmatic dispatch via exc.code."""
    content = (ROOT / "docs/api/errors.md").read_text()
    assert "exc.code" in content, (
        "errors.md missing example of programmatic dispatch via exc.code"
    )
    assert "DegenerateInput" in content, (
        "errors.md missing RegistrationError.DegenerateInput in code example"
    )


# ---------------------------------------------------------------------------
# M2: aligned_features documented in types.md
# ---------------------------------------------------------------------------


def test_types_md_has_aligned_features():
    """docs/api/types.md result types must document aligned_features field."""
    content = (ROOT / "docs/api/types.md").read_text()
    assert "aligned_features" in content, (
        "types.md missing aligned_features field documentation"
    )
