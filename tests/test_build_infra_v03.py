"""Tests for workspace-eeg.1: Build infrastructure changes for v0.3.

These tests verify the 6 changes described in bead workspace-eeg.1:
1. ci.yml: remove macos-13, update manylinux image, add test exclusion
2. CMakeLists.txt: add STABLE_ABI flag
3. pyproject.toml: version bump, docs extra, cibuildwheel config
4. meshmonk/__init__.py: importlib.metadata version
5. ADR-001: amend D2 for 12-wheel strategy
6. OpenMP verification (documentation comment in CMakeLists.txt)

TDD: these tests are written BEFORE the implementation changes.
"""

import subprocess
import sys
from pathlib import Path

WORKSPACE = Path(__file__).parent.parent

# ---------------------------------------------------------------------------
# Change 1: ci.yml changes
# ---------------------------------------------------------------------------


def test_ci_no_macos13():
    """ci.yml must not contain macos-13 runner entries."""
    ci_content = (WORKSPACE / ".github/workflows/ci.yml").read_text()
    assert (
        "macos-13" not in ci_content
    ), "ci.yml still contains macos-13 entries (deprecated December 2025)"


def test_ci_no_manylinux2014():
    """ci.yml must not reference manylinux2014 image."""
    ci_content = (WORKSPACE / ".github/workflows/ci.yml").read_text()
    assert (
        "manylinux2014" not in ci_content
    ), "ci.yml still uses manylinux2014; should be manylinux_2_28"


def test_ci_test_command_excludes_silent_logger():
    """ci.yml CIBW_TEST_COMMAND must exclude test_silent_logger.py."""
    ci_content = (WORKSPACE / ".github/workflows/ci.yml").read_text()
    assert (
        "test_silent_logger" in ci_content
    ), "ci.yml CIBW_TEST_COMMAND must add --ignore for test_silent_logger.py"


# ---------------------------------------------------------------------------
# Change 2: CMakeLists.txt STABLE_ABI flag
# ---------------------------------------------------------------------------


def test_cmake_stable_abi():
    """CMakeLists.txt must have STABLE_ABI in nanobind_add_module call."""
    cmake_content = (WORKSPACE / "CMakeLists.txt").read_text()
    assert (
        "STABLE_ABI" in cmake_content
    ), "CMakeLists.txt is missing STABLE_ABI flag in nanobind_add_module"


# ---------------------------------------------------------------------------
# Change 3: pyproject.toml changes
# ---------------------------------------------------------------------------


def test_pyproject_version_bump():
    """pyproject.toml must have version 0.3.0.dev0."""
    import tomllib

    with open(WORKSPACE / "pyproject.toml", "rb") as f:
        data = tomllib.load(f)
    assert (
        data["project"]["version"] == "0.3.0.dev0"
    ), f"Expected version 0.3.0.dev0, got {data['project']['version']}"


def test_pyproject_has_pyyaml():
    """pyproject.toml dev extra must include pyyaml."""
    import tomllib

    with open(WORKSPACE / "pyproject.toml", "rb") as f:
        data = tomllib.load(f)
    dev_deps = data["project"]["optional-dependencies"]["dev"]
    assert any(
        "pyyaml" in dep.lower() for dep in dev_deps
    ), "pyproject.toml dev extra missing pyyaml"


def test_pyproject_has_tomli():
    """pyproject.toml dev extra must include tomli (for Python < 3.11)."""
    import tomllib

    with open(WORKSPACE / "pyproject.toml", "rb") as f:
        data = tomllib.load(f)
    dev_deps = data["project"]["optional-dependencies"]["dev"]
    assert any(
        "tomli" in dep.lower() for dep in dev_deps
    ), "pyproject.toml dev extra missing tomli"


def test_pyproject_has_docs_extra():
    """pyproject.toml must have a docs optional-dependency group."""
    import tomllib

    with open(WORKSPACE / "pyproject.toml", "rb") as f:
        data = tomllib.load(f)
    assert (
        "docs" in data["project"]["optional-dependencies"]
    ), "pyproject.toml missing docs optional-dependency group"


def test_pyproject_has_cibuildwheel_section():
    """pyproject.toml must have [tool.cibuildwheel] section."""
    import tomllib

    with open(WORKSPACE / "pyproject.toml", "rb") as f:
        data = tomllib.load(f)
    assert "cibuildwheel" in data.get(
        "tool", {}
    ), "pyproject.toml missing [tool.cibuildwheel] section"


def test_pyproject_cibuildwheel_build_selector():
    """cibuildwheel build selector must include cp310, cp311, cp312."""
    import tomllib

    with open(WORKSPACE / "pyproject.toml", "rb") as f:
        data = tomllib.load(f)
    build = data["tool"]["cibuildwheel"]["build"]
    assert "cp310-*" in build
    assert "cp311-*" in build
    assert "cp312-*" in build


def test_pyproject_cibuildwheel_manylinux_image():
    """cibuildwheel linux config must use manylinux_2_28."""
    import tomllib

    with open(WORKSPACE / "pyproject.toml", "rb") as f:
        data = tomllib.load(f)
    linux = data["tool"]["cibuildwheel"]["linux"]
    assert linux["manylinux-x86_64-image"] == "manylinux_2_28"
    assert linux["manylinux-aarch64-image"] == "manylinux_2_28"


def test_pyproject_cibuildwheel_cp312_override():
    """cibuildwheel must have cp312 override for SKBUILD_WHEEL_PY_API."""
    import tomllib

    with open(WORKSPACE / "pyproject.toml", "rb") as f:
        data = tomllib.load(f)
    overrides = data["tool"]["cibuildwheel"].get("overrides", [])
    cp312_overrides = [o for o in overrides if o.get("select") == "cp312-*"]
    assert len(cp312_overrides) == 1, "Expected exactly one cp312-* override"
    env = cp312_overrides[0].get("environment", {})
    assert (
        env.get("SKBUILD_WHEEL_PY_API") == "cp312"
    ), "cp312 override missing SKBUILD_WHEEL_PY_API = cp312"


def test_pyproject_test_command_excludes_silent_logger():
    """cibuildwheel test-command must exclude test_silent_logger.py."""
    import tomllib

    with open(WORKSPACE / "pyproject.toml", "rb") as f:
        data = tomllib.load(f)
    test_cmd = data["tool"]["cibuildwheel"]["test-command"]
    assert "test_silent_logger" in test_cmd


def test_pyproject_valid_toml():
    """pyproject.toml must parse as valid TOML."""
    import tomllib

    with open(WORKSPACE / "pyproject.toml", "rb") as f:
        data = tomllib.load(f)
    # If we get here, TOML is valid
    assert data is not None


# ---------------------------------------------------------------------------
# Change 4: meshmonk/__init__.py importlib.metadata version
# ---------------------------------------------------------------------------


def test_init_uses_importlib_metadata():
    """meshmonk/__init__.py must use importlib.metadata for version."""
    init_content = (WORKSPACE / "meshmonk/__init__.py").read_text()
    assert (
        "importlib.metadata" in init_content
    ), "__init__.py still uses hardcoded __version__; should use importlib.metadata"


def test_init_no_hardcoded_version():
    """meshmonk/__init__.py must not have hardcoded version string."""
    init_content = (WORKSPACE / "meshmonk/__init__.py").read_text()
    assert (
        '__version__ = "0.2.0"' not in init_content
    ), "__init__.py still has hardcoded 0.2.0 version"


def test_init_has_fallback_version():
    """meshmonk/__init__.py must have a fallback version for uninstalled runs."""
    init_content = (WORKSPACE / "meshmonk/__init__.py").read_text()
    assert (
        "0.0.0.dev0" in init_content
    ), "__init__.py missing fallback version '0.0.0.dev0'"


def test_init_compiles():
    """meshmonk/__init__.py must have no syntax errors."""
    result = subprocess.run(
        [sys.executable, "-m", "py_compile", str(WORKSPACE / "meshmonk/__init__.py")],
        capture_output=True,
        text=True,
    )
    assert (
        result.returncode == 0
    ), f"meshmonk/__init__.py has syntax errors:\n{result.stderr}"


# ---------------------------------------------------------------------------
# Change 5: ADR-001 amendment for D2
# ---------------------------------------------------------------------------


def test_adr001_has_amendment():
    """ADR-001 must contain the v0.3 wheel strategy amendment."""
    adr_content = (
        WORKSPACE / "docs/decisions/ADR-001-meshmonk-modernization.md"
    ).read_text()
    assert (
        "Amendment (v0.3 implementation" in adr_content
    ), "ADR-001 D2 is missing the v0.3 amendment"


def test_adr001_amendment_mentions_12_wheels():
    """ADR-001 amendment must mention the 12-wheel strategy."""
    adr_content = (
        WORKSPACE / "docs/decisions/ADR-001-meshmonk-modernization.md"
    ).read_text()
    assert (
        "12 wheels" in adr_content
    ), "ADR-001 amendment must mention 12 wheels strategy"


def test_adr001_amendment_mentions_cp312_floor():
    """ADR-001 amendment must mention abi3 floor is Python 3.12."""
    adr_content = (
        WORKSPACE / "docs/decisions/ADR-001-meshmonk-modernization.md"
    ).read_text()
    assert (
        "3.12" in adr_content
    ), "ADR-001 amendment must reference Python 3.12 as abi3 floor"


# ---------------------------------------------------------------------------
# Change 6: OpenMP verification comment in CMakeLists.txt
# ---------------------------------------------------------------------------


def test_cmake_openmp_comment():
    """CMakeLists.txt must contain OpenMP verification comment."""
    cmake_content = (WORKSPACE / "CMakeLists.txt").read_text()
    assert (
        "OpenMP" in cmake_content
    ), "CMakeLists.txt missing OpenMP verification comment"
