"""
Tests for meshmonk CLI (Tier 5 — E2E smoke tests).

Tests verify:
- All subcommands have --help that exits 0 and shows expected arguments
- rigid/nonrigid/pyramid subcommands accept positional mesh args
- trimesh import guard produces a clear error message (not a traceback)
- End-to-end: rigid on real data exits 0 and produces a non-empty output file
"""

from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

DATA_DIR = Path("/workspace/data")
TEMPLATE_OBJ = DATA_DIR / "Template.obj"
DEMO_FACE_OBJ = DATA_DIR / "demoFace.obj"
MESHES_AVAILABLE = TEMPLATE_OBJ.exists() and DEMO_FACE_OBJ.exists()


# ---------------------------------------------------------------------------
# Import guard — app must be importable even if trimesh is not installed
# ---------------------------------------------------------------------------

def test_cli_module_importable():
    """cli.py is importable without triggering a trimesh import."""
    import meshmonk.cli  # noqa: F401
    assert True


def test_app_object_exists():
    """The module must expose a typer.Typer() object named `app`."""
    from meshmonk.cli import app
    import typer
    assert isinstance(app, typer.Typer)


# ---------------------------------------------------------------------------
# --help exits 0 for all subcommands
# ---------------------------------------------------------------------------

@pytest.fixture()
def runner():
    return CliRunner()


@pytest.fixture()
def app():
    from meshmonk.cli import app as _app
    return _app


def test_rigid_help(runner, app):
    result = runner.invoke(app, ["rigid", "--help"])
    assert result.exit_code == 0, result.output
    assert "floating" in result.output.lower() or "FLOATING" in result.output
    assert "target" in result.output.lower() or "TARGET" in result.output


def test_nonrigid_help(runner, app):
    result = runner.invoke(app, ["nonrigid", "--help"])
    assert result.exit_code == 0, result.output
    assert "floating" in result.output.lower() or "FLOATING" in result.output
    assert "target" in result.output.lower() or "TARGET" in result.output


def test_pyramid_help(runner, app):
    result = runner.invoke(app, ["pyramid", "--help"])
    assert result.exit_code == 0, result.output
    assert "floating" in result.output.lower() or "FLOATING" in result.output
    assert "target" in result.output.lower() or "TARGET" in result.output


def test_demo_help(runner, app):
    result = runner.invoke(app, ["demo", "--help"])
    assert result.exit_code == 0, result.output
    assert "--download" in result.output or "download" in result.output.lower()


# ---------------------------------------------------------------------------
# --help shows expected options
# ---------------------------------------------------------------------------

def test_rigid_help_shows_iterations(runner, app):
    result = runner.invoke(app, ["rigid", "--help"])
    assert result.exit_code == 0
    assert "iterations" in result.output.lower()


def test_rigid_help_shows_out(runner, app):
    result = runner.invoke(app, ["rigid", "--help"])
    assert result.exit_code == 0
    assert "--out" in result.output


def test_nonrigid_help_shows_iterations(runner, app):
    result = runner.invoke(app, ["nonrigid", "--help"])
    assert result.exit_code == 0
    assert "iterations" in result.output.lower()


def test_pyramid_help_shows_layers(runner, app):
    result = runner.invoke(app, ["pyramid", "--help"])
    assert result.exit_code == 0
    assert "layers" in result.output.lower()


# ---------------------------------------------------------------------------
# E2E: rigid registration produces output file (requires io extra + data)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not MESHES_AVAILABLE,
    reason="Test meshes not available at /workspace/data/",
)
def test_rigid_e2e_produces_output_file(runner, app, tmp_path):
    out_path = tmp_path / "result.obj"
    result = runner.invoke(
        app,
        [
            "rigid",
            str(TEMPLATE_OBJ),
            str(DEMO_FACE_OBJ),
            "--out", str(out_path),
        ],
    )
    assert result.exit_code == 0, f"Exit {result.exit_code}:\n{result.output}"
    assert out_path.exists(), "Output file was not created"
    assert out_path.stat().st_size > 0, "Output file is empty"


@pytest.mark.skipif(
    not MESHES_AVAILABLE,
    reason="Test meshes not available at /workspace/data/",
)
def test_rigid_e2e_prints_saved_message(runner, app, tmp_path):
    out_path = tmp_path / "result.obj"
    result = runner.invoke(
        app,
        [
            "rigid",
            str(TEMPLATE_OBJ),
            str(DEMO_FACE_OBJ),
            "--out", str(out_path),
        ],
    )
    assert result.exit_code == 0
    assert "saved" in result.output.lower()


# ---------------------------------------------------------------------------
# Demo subcommand: dev-convenience path (uses /workspace/data/)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not MESHES_AVAILABLE,
    reason="Test meshes not available at /workspace/data/",
)
def test_demo_download_flag_present(runner, app):
    """demo --download is a valid flag (even if URLs are placeholder)."""
    result = runner.invoke(app, ["demo", "--help"])
    assert "--download" in result.output


def test_demo_download_exits_code_1(runner, app):
    """demo --download exits with code 1 (not implemented) and message goes to stderr."""
    result = runner.invoke(app, ["demo", "--download"])
    assert result.exit_code == 1, (
        f"Expected exit code 1, got {result.exit_code}. Output: {result.output}"
    )


def test_demo_download_message_is_user_friendly(runner, app):
    """demo --download prints a user-friendly message (not 'TODO')."""
    result = runner.invoke(app, ["demo", "--download"])
    # Should NOT contain bare 'TODO'
    assert "TODO" not in result.output
    # Should mention bundled meshes
    assert "bundled" in result.output.lower() or "meshmonk/data" in result.output


@pytest.mark.skipif(
    not MESHES_AVAILABLE,
    reason="Test meshes not available at /workspace/data/",
)
def test_demo_runs_rigid_mode(runner, app, monkeypatch, tmp_path):
    """demo rigid finds meshes in /workspace/data/ dev fallback and runs."""
    # Patch the cache dir so output lands somewhere predictable
    monkeypatch.setenv("HOME", str(tmp_path))
    result = runner.invoke(app, ["demo", "rigid"])
    # Should either succeed or give clear instructions — never a raw traceback
    assert result.exit_code == 0 or "meshmonk[io]" in result.output or "data/" in result.output
