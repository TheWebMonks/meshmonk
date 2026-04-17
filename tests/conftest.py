"""
pytest conftest.py — shared fixtures and marker registration.

Custom markers:
  slow        — tests that take >30 s each (Tier 3 golden comparisons)
  integration — tests that require network access or external tools
  golden      — tests that compare against frozen golden .npz files
  advisory    — tests that report discrepancies but do not gate CI
"""

from pathlib import Path

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


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def mesh_data_dir():
    """Path to /workspace/data/ directory.

    Skips the test if the data directory is empty or does not contain any OBJ files.
    """
    data_dir = Path(__file__).parent.parent / "data"
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
