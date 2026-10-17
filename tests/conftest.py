"""
pytest conftest.py — shared fixtures and marker registration.

Custom markers:
  slow        — tests that take >30 s each (Tier 3 golden comparisons)
  integration — tests that require network access or external tools
  golden      — tests that compare against frozen golden .npz files
  advisory    — tests that report discrepancies but do not gate CI
"""

from pathlib import Path

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Marker registration
# ---------------------------------------------------------------------------


def pytest_configure(config):
    """Register custom markers to suppress PytestUnknownMarkWarning."""
    config.addinivalue_line(
        "markers",
        "slow: marks tests as slow (deselect with '-m \"not slow\"')",
    )
    config.addinivalue_line(
        "markers",
        "integration: marks tests as integration tests requiring external resources",
    )
    config.addinivalue_line(
        "markers",
        "golden: marks tests as golden-comparison tests (compare against frozen .npz files)",
    )
    config.addinivalue_line(
        "markers",
        "advisory: marks tests as advisory (failures reported but do not gate CI)",
    )


# Single repo-relative data dir; test modules import this (tracked meshes live here).
DATA_DIR = Path(__file__).parent.parent / "data"


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def mesh_data_dir():
    """Path to the repo's tracked data/ directory (DATA_DIR).

    Skips the test if the data directory is empty or does not contain any OBJ files.
    """
    data_dir = DATA_DIR
    obj_files = list(data_dir.glob("*.obj"))
    if not obj_files:
        pytest.skip(f"No OBJ files found in {data_dir} — skipping data-dependent test")
    return data_dir


@pytest.fixture
def golden_dir():
    """Path to tests/golden/ directory.

    Skips the test if the golden directory contains no .npz files.
    """
    golden = Path(__file__).parent / "golden"
    npz_files = list(golden.glob("*.npz"))
    if not npz_files:
        pytest.skip(
            f"No .npz golden files found in {golden} — run golden capture first"
        )
    return golden


@pytest.fixture
def empty_registration_kwargs():
    """rigid_register kwargs for all-empty meshes: raises MeshMonkError(DegenerateInput)."""
    empty = np.zeros((0, 6), dtype=np.float32)
    faces = np.zeros((0, 3), dtype=np.int32)
    flags = np.zeros((0,), dtype=np.float32)
    return dict(
        floating_features=empty,
        target_features=empty,
        floating_faces=faces,
        target_faces=faces,
        floating_flags=flags,
        target_flags=flags,
    )


@pytest.fixture(scope="module")
def meshes():
    """Template.obj (floating) and demoFace.obj (target) as float32/int32 arrays, once per module."""
    trimesh = pytest.importorskip("trimesh")
    template = DATA_DIR / "Template.obj"
    demoface = DATA_DIR / "demoFace.obj"
    if not template.exists() or not demoface.exists():
        pytest.skip(f"OBJ data missing in {DATA_DIR}")
    floating = trimesh.load(str(template), process=False)
    target = trimesh.load(str(demoface), process=False)
    return {
        "floating_v": np.asarray(floating.vertices, dtype=np.float32),
        "floating_f": np.asarray(floating.faces, dtype=np.int32),
        "target_v": np.asarray(target.vertices, dtype=np.float32),
        "target_f": np.asarray(target.faces, dtype=np.int32),
    }
