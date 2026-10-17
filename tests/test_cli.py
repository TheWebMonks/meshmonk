"""
Tests for meshmonk CLI (Tier 5 — E2E smoke tests).

Tests verify:
- All subcommands have --help that exits 0 and shows expected arguments
- rigid/nonrigid/pyramid subcommands accept positional mesh args
- trimesh import guard produces a clear error message (not a traceback)
- End-to-end: rigid on real data exits 0 and produces a non-empty output file
"""

from __future__ import annotations

import re
import sys

import pytest
from typer.testing import CliRunner

_ANSI_RE = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")


def _strip_ansi(s: str) -> str:
    """Strip ANSI escape sequences so substring checks survive Rich styling."""
    return _ANSI_RE.sub("", s)


from tests.conftest import DATA_DIR
TEMPLATE_OBJ = DATA_DIR / "Template.obj"
DEMOFACE_OBJ = DATA_DIR / "demoFace.obj"
MESHES_AVAILABLE = TEMPLATE_OBJ.exists() and DEMOFACE_OBJ.exists()


# ---------------------------------------------------------------------------
# Import guard — app must be importable even if trimesh is not installed
# ---------------------------------------------------------------------------


def test_cli_module_importable():
    """cli.py is importable without triggering a trimesh import."""
    import subprocess
    import sys

    code = "import sys, meshmonk.cli; sys.exit('trimesh' in sys.modules)"
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr or "importing meshmonk.cli imported trimesh"


def test_require_trimesh_missing_gives_clear_error(monkeypatch, capsys):
    """Without trimesh, _require_trimesh exits 1 with the install hint."""
    import click

    import meshmonk.cli as cli

    monkeypatch.setitem(sys.modules, "trimesh", None)  # makes `import trimesh` raise
    with pytest.raises((SystemExit, click.exceptions.Exit)) as ei:
        cli._require_trimesh()
    assert ei.value.exit_code == 1
    assert "trimesh is required" in capsys.readouterr().err


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
    # option row (flag followed by its help column), not the docstring prose
    assert re.search(r"^\W*--download\s{2,}", _strip_ansi(result.output), re.M)


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
    assert "--out" in _strip_ansi(result.output)


def test_nonrigid_help_shows_iterations(runner, app):
    result = runner.invoke(app, ["nonrigid", "--help"])
    assert result.exit_code == 0
    assert "iterations" in result.output.lower()


def test_pyramid_help_shows_layers(runner, app):
    result = runner.invoke(app, ["pyramid", "--help"])
    assert result.exit_code == 0
    assert re.search(r"^\W*--layers\s{2,}", _strip_ansi(result.output), re.M)


# ---------------------------------------------------------------------------
# E2E: rigid registration produces output file (requires io extra + data)
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not MESHES_AVAILABLE,
    reason="Test meshes not available at data/",
)
def test_rigid_e2e_produces_output_file(runner, app, tmp_path):
    out_path = tmp_path / "result.obj"
    result = runner.invoke(
        app,
        [
            "rigid",
            str(TEMPLATE_OBJ),
            str(DEMOFACE_OBJ),
            "--out",
            str(out_path),
        ],
    )
    assert result.exit_code == 0, f"Exit {result.exit_code}:\n{result.output}"
    assert out_path.exists(), "Output file was not created"
    assert out_path.stat().st_size > 0, "Output file is empty"


@pytest.mark.skipif(
    not MESHES_AVAILABLE,
    reason="Test meshes not available at data/",
)
def test_rigid_e2e_prints_saved_message(runner, app, tmp_path):
    out_path = tmp_path / "result.obj"
    result = runner.invoke(
        app,
        [
            "rigid",
            str(TEMPLATE_OBJ),
            str(DEMOFACE_OBJ),
            "--out",
            str(out_path),
        ],
    )
    assert result.exit_code == 0
    assert "saved" in result.output.lower()


# ---------------------------------------------------------------------------
# Demo subcommand: dev-convenience path (uses data/)
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not MESHES_AVAILABLE,
    reason="Test meshes not available at data/",
)
def test_demo_download_flag_present(runner, app):
    """demo --download is a valid flag (even if URLs are placeholder)."""
    result = runner.invoke(app, ["demo", "--help"])
    assert "--download" in result.output


def test_demo_download_exits_code_1(runner, app):
    """demo --download exits with code 1 in environments where GitHub Releases are unreachable."""
    result = runner.invoke(app, ["demo", "--download"])
    assert (
        result.exit_code == 1
    ), f"Expected exit code 1, got {result.exit_code}. Output: {result.output}"


def test_demo_download_message_is_user_friendly(runner, app):
    """demo --download prints a user-friendly message (not 'TODO')."""
    result = runner.invoke(app, ["demo", "--download"])
    # Should NOT contain bare 'TODO'
    assert "TODO" not in result.output
    # Should attempt to download (new implementation uses real GitHub Releases URLs)
    assert "Downloading" in result.output or "download" in result.output.lower()


@pytest.mark.skipif(
    not MESHES_AVAILABLE,
    reason="Test meshes not available at data/",
)
def test_demo_runs_rigid_mode(runner, app, monkeypatch, tmp_path):
    """demo rigid finds meshes via the repo data/ fallback, runs, and writes output."""
    import meshmonk.cli as cli

    # _CACHE_DIR is bound at import; patch it so the empty cache falls through
    # to cli._REPO_DATA_DIR and demo_result.obj lands in tmp_path
    monkeypatch.setattr(cli, "_CACHE_DIR", tmp_path)
    result = runner.invoke(app, ["demo", "rigid"])
    assert result.exit_code == 0, result.output
    assert (tmp_path / "demo_result.obj").stat().st_size > 0
