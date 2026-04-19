"""
Tier 5 E2E smoke tests.

Verifies that the CLI entry point is installed and functional.
These tests are deliberately minimal — they check exit codes and output existence,
not numerical correctness (which is covered by Tiers 1-4).
"""

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

DATA_DIR = Path(__file__).parent.parent / "data"
TEMPLATE_OBJ = DATA_DIR / "Template.obj"
DEMOFACE_OBJ = DATA_DIR / "demoFace.obj"


def _meshmonk_cmd() -> list:
    """Return the command prefix for invoking the meshmonk CLI.

    Prefers the installed 'meshmonk' entry point if found on PATH.
    Falls back to 'python -m meshmonk.cli' using the current interpreter.
    """
    # Look for the installed entry point next to the current Python executable
    # (i.e. in the same venv bin/ directory) before falling back to PATH search.
    venv_bin = Path(sys.executable).parent
    venv_meshmonk = venv_bin / "meshmonk"
    if venv_meshmonk.exists():
        return [str(venv_meshmonk)]

    # Fall back to PATH-based search
    cli = shutil.which("meshmonk")
    if cli:
        return [cli]

    # Last resort: use typer's __main__ invocation
    return [sys.executable, "-m", "meshmonk.cli"]


def test_cli_help_exits_zero():
    """Assert 'meshmonk --help' exits with code 0."""
    cmd = _meshmonk_cmd() + ["--help"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, (
        f"'meshmonk --help' exited with code {result.returncode}.\n"
        f"stdout: {result.stdout[:500]}\n"
        f"stderr: {result.stderr[:500]}"
    )


def test_cli_rigid_exits_zero(tmp_path):
    """Assert 'meshmonk rigid' completes and produces a non-empty output file.

    Skipped if:
      - trimesh is not installed (required for OBJ I/O by the CLI)
      - Template.obj or demoFace.obj are not present in data/
    """
    pytest.importorskip("trimesh")

    if not TEMPLATE_OBJ.exists():
        pytest.skip(f"Template.obj not found at {TEMPLATE_OBJ}")
    if not DEMOFACE_OBJ.exists():
        pytest.skip(f"demoFace.obj not found at {DEMOFACE_OBJ}")

    out_file = tmp_path / "result.obj"
    cmd = _meshmonk_cmd() + [
        "rigid",
        str(TEMPLATE_OBJ),
        str(DEMOFACE_OBJ),
        "--out",
        str(out_file),
    ]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=300,  # rigid registration can take up to 180 s in CI
    )
    assert result.returncode == 0, (
        f"'meshmonk rigid' exited with code {result.returncode}.\n"
        f"stdout: {result.stdout[:500]}\n"
        f"stderr: {result.stderr[:500]}"
    )
    assert out_file.exists(), f"Output file {out_file} was not created"
    assert out_file.stat().st_size > 0, f"Output file {out_file} is empty"
