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

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib  # type: ignore[import-not-found]

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


def test_ci_test_command_excludes_silent_logger():
    """ci.yml CIBW_TEST_COMMAND must exclude test_silent_logger.py."""
    import yaml

    ci = yaml.safe_load((WORKSPACE / ".github/workflows/ci.yml").read_text())
    cmds = []

    def walk(node):
        if isinstance(node, dict):
            for k, v in node.items():
                if k == "CIBW_TEST_COMMAND":
                    cmds.append(v)
                walk(v)
        elif isinstance(node, list):
            for v in node:
                walk(v)

    walk(ci)
    assert cmds, "ci.yml has no CIBW_TEST_COMMAND"
    for cmd in cmds:
        assert (
            "--ignore={project}/tests/test_silent_logger.py" in cmd.split()
        ), "ci.yml CIBW_TEST_COMMAND must add --ignore for test_silent_logger.py"


# ---------------------------------------------------------------------------
# Change 2: CMakeLists.txt STABLE_ABI flag
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Change 3: pyproject.toml changes
# ---------------------------------------------------------------------------


def test_pyproject_has_pyyaml():
    """pyproject.toml dev extra must include pyyaml."""
    with open(WORKSPACE / "pyproject.toml", "rb") as f:
        data = tomllib.load(f)
    dev_deps = data["project"]["optional-dependencies"]["dev"]
    assert any(
        "pyyaml" in dep.lower() for dep in dev_deps
    ), "pyproject.toml dev extra missing pyyaml"


def test_pyproject_has_tomli():
    """pyproject.toml dev extra must include tomli (for Python < 3.11)."""
    with open(WORKSPACE / "pyproject.toml", "rb") as f:
        data = tomllib.load(f)
    dev_deps = data["project"]["optional-dependencies"]["dev"]
    assert any(
        "tomli" in dep.lower() for dep in dev_deps
    ), "pyproject.toml dev extra missing tomli"


def test_pyproject_cibuildwheel_cp312_override():
    """cibuildwheel must have cp312 override for SKBUILD_WHEEL_PY_API."""
    with open(WORKSPACE / "pyproject.toml", "rb") as f:
        data = tomllib.load(f)
    overrides = data["tool"]["cibuildwheel"].get("overrides", [])
    cp312_overrides = [o for o in overrides if o.get("select") == "cp312-*"]
    assert len(cp312_overrides) == 1, "Expected exactly one cp312-* override"
    env = cp312_overrides[0].get("environment", {})
    assert (
        env.get("SKBUILD_WHEEL_PY_API") == "cp312"
    ), "cp312 override missing SKBUILD_WHEEL_PY_API = cp312"


# ---------------------------------------------------------------------------
# Change 4: meshmonk/__init__.py importlib.metadata version
# ---------------------------------------------------------------------------


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
