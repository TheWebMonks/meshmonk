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


class _FakeResponse:
    """Stand-in for the urlopen context manager: serves fixed bytes."""

    def __init__(self, data: bytes):
        self._data = data

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def read(self) -> bytes:
        return self._data


@pytest.fixture
def offline_cache(monkeypatch, tmp_path):
    """Point the demo cache at tmp_path and make every urlopen fail offline."""
    import urllib.error
    import urllib.request

    import meshmonk.cli as cli

    def _unreachable(*args, **kwargs):
        raise urllib.error.URLError("network disabled in tests")

    monkeypatch.setattr(cli, "_CACHE_DIR", tmp_path)
    monkeypatch.setattr(urllib.request, "urlopen", _unreachable)
    return tmp_path


def _serve(monkeypatch, body_for):
    """Patch urlopen to return body_for(name) for each DEMO_ASSETS url."""
    import urllib.request

    import meshmonk.cli as cli

    by_url = {a["url"]: name for name, a in cli.DEMO_ASSETS.items()}
    monkeypatch.setattr(
        urllib.request,
        "urlopen",
        lambda url, timeout=None: _FakeResponse(body_for(by_url[url])),
    )


def test_demo_download_exits_code_1(runner, app, offline_cache):
    """demo --download exits with code 1 when GitHub Releases are unreachable."""
    result = runner.invoke(app, ["demo", "--download"])
    assert (
        result.exit_code == 1
    ), f"Expected exit code 1, got {result.exit_code}. Output: {result.output}"


def test_demo_download_message_is_user_friendly(runner, app, offline_cache):
    """demo --download prints a user-friendly message (not 'TODO')."""
    result = runner.invoke(app, ["demo", "--download"])
    # Should NOT contain bare 'TODO'
    assert "TODO" not in result.output
    assert "Downloading" in result.output or "download" in result.output.lower()


@pytest.mark.parametrize("name", ["Template.obj", "demoFace.obj"])
def test_demo_asset_digest_matches_tracked_file(name):
    """DEMO_ASSETS pins the sha256 of the tracked data/ file it ships as."""
    import hashlib

    import meshmonk.cli as cli

    expected = hashlib.sha256((DATA_DIR / name).read_bytes()).hexdigest()
    assert cli.DEMO_ASSETS[name]["sha256"] == expected


def test_demo_download_digest_mismatch_exits_1_and_removes_file(
    runner, app, monkeypatch, tmp_path
):
    """Wrong bytes from the release URL: exit 1 and no file left in the cache."""
    import meshmonk.cli as cli

    monkeypatch.setattr(cli, "_CACHE_DIR", tmp_path)
    _serve(monkeypatch, lambda name: b"not the real mesh")
    result = runner.invoke(app, ["demo", "--download"])
    assert result.exit_code == 1, result.output
    assert "SHA-256 mismatch for Template.obj" in result.output
    assert not (tmp_path / "Template.obj").exists()


def test_demo_download_writes_verified_assets(runner, app, monkeypatch, tmp_path):
    """Correct bytes from the release URLs: exit 0, both files cached intact."""
    import meshmonk.cli as cli

    monkeypatch.setattr(cli, "_CACHE_DIR", tmp_path)
    _serve(monkeypatch, lambda name: (DATA_DIR / name).read_bytes())
    result = runner.invoke(app, ["demo", "--download"])
    assert result.exit_code == 0, result.output
    assert "integrity check not configured" not in result.output
    for name in cli.DEMO_ASSETS:
        assert (tmp_path / name).read_bytes() == (DATA_DIR / name).read_bytes()


def test_demo_not_found_hint_names_real_lookup_dirs(runner, app, monkeypatch):
    """The not-found hint names the cache and package data/ dirs, not the cwd (c75)."""
    import meshmonk.cli as cli

    monkeypatch.setattr(cli, "_find_demo_meshes", lambda: None)
    result = runner.invoke(app, ["demo", "rigid"])
    out = _strip_ansi(result.output)
    assert result.exit_code == 1, out
    assert "~/.cache/meshmonk/" in out
    assert "data/ directory next to the meshmonk package" in out
    assert "repository root" not in out


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
