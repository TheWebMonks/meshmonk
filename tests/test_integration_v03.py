"""Integration tests for v0.3 configuration artifacts.

Verifies pyproject.toml, ci.yml, release.yml, mkdocs.yml, and source files are
correctly configured for the v0.3 PyPI-ready release. No live builds or network
access required.

Focus: CROSS-CUTTING concerns — verifying that the pieces work TOGETHER, not
re-testing individual file changes already covered by test_build_infra_v03.py.

Requires: pyyaml (pip install pyyaml or pip install ".[dev]")
tomllib is stdlib in Python 3.11+; on Python 3.10 install tomli:
  pip install "tomli>=2.0"  or  pip install ".[dev]" (tomli is in dev extras)
"""
from __future__ import annotations

import re
from pathlib import Path

import yaml  # pip install pyyaml

try:
    import tomllib
except ImportError:
    import tomli as tomllib  # type: ignore[no-redef]  # Python 3.10 fallback

ROOT = Path(__file__).parent.parent


# ---------------------------------------------------------------------------
# Section 1: Build config cross-cutting correctness
# (release workflow references the same cibuildwheel config in pyproject.toml)
# ---------------------------------------------------------------------------


def test_pyproject_toml_is_valid_toml():
    """pyproject.toml must parse without error."""
    content = (ROOT / "pyproject.toml").read_text()
    tomllib.loads(content)  # raises on invalid TOML


def test_pyproject_version_is_v03():
    """Version must be 0.3.x (dev, rc, or final)."""
    content = tomllib.loads((ROOT / "pyproject.toml").read_text())
    version = content["project"]["version"]
    assert version.startswith("0.3."), f"Expected 0.3.x, got {version!r}"


def test_pyproject_has_cibuildwheel_config():
    """pyproject.toml must have [tool.cibuildwheel] section."""
    content = tomllib.loads((ROOT / "pyproject.toml").read_text())
    assert "cibuildwheel" in content.get("tool", {}), "Missing [tool.cibuildwheel]"


def test_cibuildwheel_build_selectors():
    """Build selector must include cp310, cp311, cp312."""
    content = tomllib.loads((ROOT / "pyproject.toml").read_text())
    build = content["tool"]["cibuildwheel"]["build"]
    assert "cp310-*" in build
    assert "cp311-*" in build
    assert "cp312-*" in build


def test_cibuildwheel_manylinux_image_is_2_28():
    """manylinux image must be manylinux_2_28, not manylinux2014."""
    content = tomllib.loads((ROOT / "pyproject.toml").read_text())
    linux_cfg = content["tool"]["cibuildwheel"]["linux"]
    assert linux_cfg["manylinux-x86_64-image"] == "manylinux_2_28"
    assert linux_cfg["manylinux-aarch64-image"] == "manylinux_2_28"


def test_cibuildwheel_macos_arm64_only():
    """macOS must target arm64 only (no x86_64 — runner deprecated December 2025)."""
    content = tomllib.loads((ROOT / "pyproject.toml").read_text())
    macos_archs = content["tool"]["cibuildwheel"]["macos"]["archs"]
    assert macos_archs == ["arm64"]
    assert "x86_64" not in macos_archs


def test_cibuildwheel_cp312_override_sets_abi3():
    """cp312 override must set SKBUILD_WHEEL_PY_API=cp312 for abi3 tagging."""
    content = tomllib.loads((ROOT / "pyproject.toml").read_text())
    overrides = content["tool"]["cibuildwheel"]["overrides"]
    cp312_override = next((o for o in overrides if o["select"] == "cp312-*"), None)
    assert cp312_override is not None, "Missing cp312-* override"
    env = cp312_override.get("environment", {})
    assert env.get("SKBUILD_WHEEL_PY_API") == "cp312"


def test_cibuildwheel_test_excludes_silent_logger():
    """cibuildwheel test-command must exclude test_silent_logger.py."""
    content = tomllib.loads((ROOT / "pyproject.toml").read_text())
    test_cmd = content["tool"]["cibuildwheel"]["test-command"]
    assert "test_silent_logger" in test_cmd


def test_no_py_api_in_scikit_build_wheel():
    """Must NOT have py-api in [tool.scikit-build.wheel] (would tag all wheels as abi3)."""
    content = tomllib.loads((ROOT / "pyproject.toml").read_text())
    skbuild = content.get("tool", {}).get("scikit-build", {})
    wheel = skbuild.get("wheel", {})
    assert "py-api" not in wheel, "py-api in scikit-build.wheel would break cp310/cp311 tagging"


def test_docs_extra_in_pyproject():
    """pyproject.toml must have docs optional-dependency group."""
    content = tomllib.loads((ROOT / "pyproject.toml").read_text())
    extras = content["project"]["optional-dependencies"]
    assert "docs" in extras
    docs = extras["docs"]
    assert any("mkdocs-material" in d for d in docs)
    assert any("mkdocstrings" in d for d in docs)


# ---------------------------------------------------------------------------
# Section 2: CI workflow + pyproject.toml consistency
# (ci.yml dry-run must use same manylinux image as pyproject.toml release config)
# ---------------------------------------------------------------------------


def test_ci_yml_has_no_macos13():
    """ci.yml must not reference macos-13 (runner deprecated December 2025)."""
    content = (ROOT / ".github/workflows/ci.yml").read_text()
    assert "macos-13" not in content, "macos-13 runner still referenced in ci.yml"


def test_ci_yml_and_pyproject_agree_on_manylinux_image():
    """ci.yml and pyproject.toml must reference the same manylinux image version."""
    ci_content = (ROOT / ".github/workflows/ci.yml").read_text()
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text())
    pyproject_image = pyproject["tool"]["cibuildwheel"]["linux"]["manylinux-x86_64-image"]
    # ci.yml must reference the same image used in the release config
    assert pyproject_image in ci_content, (
        f"ci.yml does not reference pyproject.toml manylinux image '{pyproject_image}'"
    )


def test_ci_yml_no_manylinux2014_while_pyproject_uses_2_28():
    """ci.yml must not use manylinux2014 when pyproject.toml specifies manylinux_2_28."""
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text())
    manylinux_image = pyproject["tool"]["cibuildwheel"]["linux"]["manylinux-x86_64-image"]
    assert manylinux_image == "manylinux_2_28", "pyproject.toml image is not manylinux_2_28"
    ci_content = (ROOT / ".github/workflows/ci.yml").read_text()
    assert "manylinux2014" not in ci_content, (
        "ci.yml still uses manylinux2014 — inconsistent with pyproject.toml"
    )


# ---------------------------------------------------------------------------
# Section 3: Release workflow correctness and cross-file consistency
# (release.yml uses cibuildwheel from pyproject.toml; environments exist)
# ---------------------------------------------------------------------------


def test_release_yml_exists():
    """release.yml must exist."""
    assert (ROOT / ".github/workflows/release.yml").exists()


def test_release_yml_is_valid_yaml():
    """release.yml must parse as valid YAML."""
    content = (ROOT / ".github/workflows/release.yml").read_text()
    yaml.safe_load(content)  # raises on invalid YAML


def test_release_yml_has_trusted_publisher_permissions():
    """release.yml must have id-token: write and contents: write."""
    content = (ROOT / ".github/workflows/release.yml").read_text()
    assert "id-token: write" in content
    assert "contents: write" in content


def test_release_yml_has_pypi_environment():
    """release.yml must reference pypi-release environment."""
    content = (ROOT / ".github/workflows/release.yml").read_text()
    assert "pypi-release" in content


def test_release_yml_has_testpypi_environment():
    """release.yml must reference testpypi-release environment."""
    content = (ROOT / ".github/workflows/release.yml").read_text()
    assert "testpypi-release" in content


def test_release_yml_has_workflow_dispatch():
    """release.yml must have workflow_dispatch trigger."""
    content = (ROOT / ".github/workflows/release.yml").read_text()
    parsed = yaml.safe_load(content)
    # YAML parses bare 'on' as Python True (truthy value), so check both keys
    triggers = parsed.get("on", parsed.get(True, {}))
    assert "workflow_dispatch" in triggers


def test_release_yml_actions_are_sha_pinned():
    """All 'uses:' lines in release.yml must reference a 40-char lowercase hex SHA."""
    content = (ROOT / ".github/workflows/release.yml").read_text()
    uses_lines = re.findall(r"uses:\s+(\S+)", content)
    for ref in uses_lines:
        sha_part = ref.split("@")[-1] if "@" in ref else ""
        assert (
            len(sha_part) == 40
            and re.fullmatch(r"[0-9a-f]{40}", sha_part) is not None
        ), f"Action not SHA-pinned with 40-char lowercase hex: {ref}"


def test_release_yml_references_cibuildwheel():
    """release.yml must use pypa/cibuildwheel to build wheels."""
    content = (ROOT / ".github/workflows/release.yml").read_text()
    assert "cibuildwheel" in content, "release.yml does not reference cibuildwheel"


def test_release_yml_uses_pypi_publish_action():
    """release.yml must use pypa/gh-action-pypi-publish for uploads."""
    content = (ROOT / ".github/workflows/release.yml").read_text()
    assert "pypa/gh-action-pypi-publish" in content


# ---------------------------------------------------------------------------
# Section 4: Version management — single source of truth
# (__init__.py reads from pyproject.toml via importlib.metadata)
# ---------------------------------------------------------------------------


def test_init_py_uses_importlib_metadata():
    """__init__.py must use importlib.metadata for __version__, not a hardcoded string."""
    content = (ROOT / "meshmonk/__init__.py").read_text()
    assert "importlib.metadata" in content
    assert 'version("meshmonk")' in content or "version('meshmonk')" in content


def test_init_py_has_fallback_version():
    """__init__.py must have a fallback version for uninstalled dev runs."""
    content = (ROOT / "meshmonk/__init__.py").read_text()
    assert "0.0.0.dev0" in content


def test_init_py_no_hardcoded_version():
    """__init__.py must not have __version__ = '0.2.0' or similar hardcoded release version."""
    content = (ROOT / "meshmonk/__init__.py").read_text()
    # The fallback "0.0.0.dev0" is intentionally allowed; match only clean X.Y.Z without suffix
    hardcoded = re.search(r'__version__\s*=\s*["\']0\.\d+\.\d+["\']', content)
    assert hardcoded is None, f"Hardcoded version found: {hardcoded.group()}"


def test_pyproject_version_is_consistent_with_v03_branch():
    """pyproject.toml version must be 0.3.x, consistent with __init__.py fallback being 0.0.0.dev0."""
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text())
    version = pyproject["project"]["version"]
    # Must be 0.3.x (not stuck at 0.2.x)
    assert version.startswith("0.3."), f"pyproject.toml still at old version: {version}"
    # __init__.py fallback must NOT be the release version (that would be dual-maintenance)
    init_content = (ROOT / "meshmonk/__init__.py").read_text()
    assert version not in init_content or "importlib.metadata" in init_content, (
        "__init__.py hardcodes the pyproject version instead of using importlib.metadata"
    )


# ---------------------------------------------------------------------------
# Section 5: Docs cross-cutting — mkdocs.yml nav references existing files
# ---------------------------------------------------------------------------


def test_mkdocs_yml_exists():
    """mkdocs.yml must exist at repo root."""
    assert (ROOT / "mkdocs.yml").exists()


def test_mkdocs_yml_is_valid_yaml():
    """mkdocs.yml must parse as valid YAML."""
    content = (ROOT / "mkdocs.yml").read_text()
    yaml.safe_load(content)


def test_mkdocs_uses_material_theme():
    """mkdocs.yml must use mkdocs-material theme."""
    content = yaml.safe_load((ROOT / "mkdocs.yml").read_text())
    assert content.get("theme", {}).get("name") == "material"


def test_docs_index_exists():
    """docs/index.md must exist."""
    assert (ROOT / "docs/index.md").exists()


def test_docs_quickstart_exists():
    """docs/quickstart.md must exist."""
    assert (ROOT / "docs/quickstart.md").exists()


def test_mkdocs_nav_references_existing_files():
    """Every file referenced in mkdocs.yml nav must exist under docs/."""
    mkdocs = yaml.safe_load((ROOT / "mkdocs.yml").read_text())
    nav = mkdocs.get("nav", [])

    def collect_paths(nav_item):
        """Recursively collect file paths from a nav entry."""
        paths = []
        if isinstance(nav_item, dict):
            for value in nav_item.values():
                paths.extend(collect_paths(value))
        elif isinstance(nav_item, list):
            for item in nav_item:
                paths.extend(collect_paths(item))
        elif isinstance(nav_item, str):
            # It's a file path (not a URL)
            if not nav_item.startswith("http"):
                paths.append(nav_item)
        return paths

    missing = []
    for path_str in collect_paths(nav):
        full_path = ROOT / "docs" / path_str
        if not full_path.exists():
            missing.append(path_str)

    assert not missing, (
        f"mkdocs.yml nav references files that don't exist under docs/: {missing}"
    )


def test_mkdocs_nav_has_decisions_section():
    """mkdocs.yml nav must reference ADR decisions doc."""
    mkdocs = yaml.safe_load((ROOT / "mkdocs.yml").read_text())
    content_str = str(mkdocs)
    assert "decisions" in content_str.lower() or "ADR" in content_str, (
        "mkdocs.yml nav does not reference decisions/ADR docs"
    )


# ---------------------------------------------------------------------------
# Section 6: CLI demo --download cross-cutting
# (DEMO_ASSETS uses /releases/latest/download/ consistent with GitHub Releases)
# ---------------------------------------------------------------------------


def test_cli_demo_assets_dict_exists():
    """cli.py must define DEMO_ASSETS dict with GitHub Releases URLs."""
    content = (ROOT / "meshmonk/cli.py").read_text()
    assert "DEMO_ASSETS" in content
    assert "releases/latest/download" in content


def test_cli_demo_assets_have_sha256_field():
    """Each DEMO_ASSETS entry must have a sha256 field (may be None)."""
    content = (ROOT / "meshmonk/cli.py").read_text()
    assert "sha256" in content


def test_cli_no_none_url_constants():
    """Old _TEMPLATE_URL = None and _DEMO_FACE_URL = None pattern must be removed."""
    content = (ROOT / "meshmonk/cli.py").read_text()
    assert "_TEMPLATE_URL" not in content
    assert "_DEMO_FACE_URL" not in content


def test_cli_demo_assets_urls_point_to_correct_repo():
    """DEMO_ASSETS URLs must reference the jsnyde0/meshmonk GitHub repo."""
    content = (ROOT / "meshmonk/cli.py").read_text()
    assert "github.com/jsnyde0/meshmonk" in content, (
        "DEMO_ASSETS URLs do not reference the expected GitHub repo"
    )


# ---------------------------------------------------------------------------
# Section 7: STABLE_ABI in CMakeLists.txt (cross-cut with pyproject.toml cp312 override)
# ---------------------------------------------------------------------------


def test_cmake_stable_abi_flag():
    """CMakeLists.txt nanobind_add_module call must include STABLE_ABI flag."""
    content = (ROOT / "CMakeLists.txt").read_text()
    assert "STABLE_ABI" in content
    match = re.search(r"nanobind_add_module\([^)]*STABLE_ABI[^)]*\)", content)
    assert match is not None, "STABLE_ABI not found inside nanobind_add_module(...)"


def test_cmake_stable_abi_consistent_with_cp312_override():
    """CMakeLists.txt STABLE_ABI flag must co-exist with cp312 override in pyproject.toml.

    If CMakeLists.txt has STABLE_ABI, pyproject.toml must have the cp312 override
    (to avoid abi3 tagging cp310/cp311 wheels). If the override is missing, the
    STABLE_ABI flag on its own would still work (nanobind ignores it on 3.10/3.11),
    but the wheels won't be tagged abi3 for 3.12+.
    """
    cmake_content = (ROOT / "CMakeLists.txt").read_text()
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text())

    has_stable_abi = "STABLE_ABI" in cmake_content
    overrides = pyproject.get("tool", {}).get("cibuildwheel", {}).get("overrides", [])
    has_cp312_override = any(o.get("select") == "cp312-*" for o in overrides)

    if has_stable_abi:
        assert has_cp312_override, (
            "CMakeLists.txt has STABLE_ABI but pyproject.toml is missing the cp312-* "
            "override (SKBUILD_WHEEL_PY_API=cp312). abi3 wheels won't be produced for 3.12+."
        )


# ---------------------------------------------------------------------------
# Section 8: displacement_field docstring in bindings.cpp
# ---------------------------------------------------------------------------


def test_displacement_field_has_docstring_in_bindings():
    """bindings.cpp displacement_field binding must have a docstring."""
    content = (ROOT / "bindings/bindings.cpp").read_text()
    # Find displacement_field def and verify it has a string argument (docstring)
    match = re.search(
        r'\.def_ro\(\s*"displacement_field"[^;]*"[^"]+"\s*\)',
        content,
        re.DOTALL,
    )
    assert match is not None, (
        "displacement_field binding in bindings.cpp is missing a docstring argument"
    )


# ---------------------------------------------------------------------------
# Section 9: ADR-001 amendment present
# ---------------------------------------------------------------------------


def test_adr001_d2_amendment_present():
    """ADR-001 must contain the D2 amendment for the v0.3 12-wheel strategy."""
    adr_path = ROOT / "docs/decisions/ADR-001-meshmonk-modernization.md"
    assert adr_path.exists(), f"ADR-001 not found at {adr_path}"
    content = adr_path.read_text()
    assert "Amendment" in content, "ADR-001 D2 amendment is missing"
    assert "12 wheels" in content, "ADR-001 amendment must mention 12-wheel strategy"
    assert "3.12" in content, "ADR-001 amendment must reference Python 3.12 as abi3 floor"
