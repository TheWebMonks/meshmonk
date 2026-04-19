"""
Integration tests for v0.2 features. Each test exercises multiple beads together.
Run after all v0.2 beads are complete.
"""
import pathlib
import subprocess
import sys

import numpy as np
import pytest

import meshmonk

# ---- Bead 1: MeshMonkError is a typed exception ----

def test_meshmonk_error_is_subclass_of_runtime_error():
    assert issubclass(meshmonk.MeshMonkError, RuntimeError)

def test_meshmonk_error_not_raised_for_valid_input():
    pytest.importorskip("trimesh")
    import trimesh
    m = trimesh.creation.icosphere(subdivisions=0)
    result = meshmonk.rigid_register(floating=m, target=m, num_iterations=2)
    assert result is not None

def test_degenerate_input_raises_meshmonk_error_not_generic_runtime_error():
    # After bead 1: MeshMonkError is only thrown for intentional library failures
    empty = np.zeros((0, 6), dtype=np.float32)
    faces = np.zeros((0, 3), dtype=np.int32)
    flags = np.zeros((0,), dtype=np.float32)
    with pytest.raises(meshmonk.MeshMonkError) as exc_info:
        meshmonk.rigid_register(
            floating_features=empty, target_features=empty,
            floating_faces=faces, target_faces=faces,
            floating_flags=flags, target_flags=flags,
        )
    # Error message must identify the error type
    assert "DegenerateInput" in str(exc_info.value)

# ---- Bead 2: set_log_level("silent") suppresses all output ----
#
# IMPORTANT: meshmonk::log() writes to C-level stderr via std::cerr.
# Python's capsys fixture does NOT capture C-level stderr (it only captures
# Python-level sys.stderr writes). Use subprocess to capture C-level output.

def test_silent_log_level_produces_no_output():
    """set_log_level('silent') must suppress ALL output including from internal classes."""
    script = (
        "import meshmonk; import trimesh; "
        "meshmonk.set_log_level('silent'); "
        "m = trimesh.creation.icosphere(subdivisions=0); "
        "meshmonk.rigid_register(floating=m, target=m, num_iterations=2)"
    )
    result = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"Subprocess failed: {result.stderr}"
    assert result.stderr == "", (
        "Expected no stderr output in silent mode, but got:\n" + result.stderr
    )

# ---- Bead 3: rigid_params kwarg works end-to-end ----

def test_pyramid_register_with_rigid_params():
    pytest.importorskip("trimesh")
    import trimesh
    m = trimesh.creation.icosphere(subdivisions=0)
    result = meshmonk.pyramid_register(
        floating=m,
        target=m,
        rigid_params={},
        num_iterations=4,
        num_pyramid_layers=1,
    )
    assert result.aligned_features.shape[1] == 6
    assert result.aligned_features.shape[0] == len(m.vertices)

def test_nonrigid_register_with_rigid_params():
    pytest.importorskip("trimesh")
    import trimesh
    m = trimesh.creation.icosphere(subdivisions=0)
    result = meshmonk.nonrigid_register(
        floating=m,
        target=m,
        rigid_params={},
        num_iterations=3,
    )
    assert result.aligned_features.shape[1] == 6

# ---- rigid_params TypeError validation ----

def test_rigid_params_type_error():
    pytest.importorskip("trimesh")
    import trimesh
    m = trimesh.creation.icosphere(subdivisions=0)
    with pytest.raises(TypeError, match="rigid_params must be a dict"):
        meshmonk.nonrigid_register(floating=m, target=m, rigid_params=True)

# ---- Bead 4: CLI demo --download exits 1 ----

def test_demo_download_exits_nonzero():
    proc = subprocess.run(
        [sys.executable, "-m", "meshmonk", "demo", "--download"],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 1, (
        f"Expected exit code 1 for 'demo --download', got {proc.returncode}. "
        f"stdout: {proc.stdout!r}, stderr: {proc.stderr!r}"
    )
    combined = proc.stdout + proc.stderr
    assert "Download URLs not yet configured" in combined, (
        "Expected message 'Download URLs not yet configured' in output. "
        f"Got stdout: {proc.stdout!r}, stderr: {proc.stderr!r}"
    )

# ---- Bead 5: py.typed marker exists ----

def test_py_typed_marker_exists():
    package_dir = pathlib.Path(meshmonk.__file__).parent
    py_typed = package_dir / "py.typed"
    assert py_typed.exists(), (
        f"py.typed missing at {py_typed}. "
        "Create an empty file at /workspace/meshmonk/py.typed"
    )

# ---- Version check ----

def test_version_is_string():
    # Version is now read from importlib.metadata (pyproject.toml is the single source
    # of truth as of v0.3). The fallback for uninstalled dev runs is "0.0.0.dev0".
    assert isinstance(meshmonk.__version__, str), (
        f"__version__ must be a string, got {type(meshmonk.__version__)}"
    )
    assert meshmonk.__version__, "__version__ must not be empty"
    assert meshmonk.__version__ != "0.2.0", (
        "__version__ is still hardcoded 0.2.0; should use importlib.metadata"
    )
