"""
Tier 4 soft anchor tests.

Verifies the existence and parseability of committed reference files in data/.
These tests act as a soft anchor — they should always pass since the files are
committed and do not require any build or registration to run.
"""

from pathlib import Path

import numpy as np
import pytest

DATA_DIR = Path(__file__).parent.parent / "data"


def test_rigid_transform_file_exists():
    """Assert data/rigid_transform.txt exists and is parseable as a 4x4 float matrix.

    This file was captured from a meshmonk_cli rigid_reg run and committed in v0.0.
    It serves as a soft numerical anchor (loose-tolerance sanity reference) for Tier 3.5.
    """
    rigid_transform_path = DATA_DIR / "rigid_transform.txt"

    assert rigid_transform_path.exists(), (
        f"data/rigid_transform.txt not found at {rigid_transform_path}. "
        "This file should be committed to the repo as part of the v0.0 legacy data."
    )

    matrix = np.loadtxt(str(rigid_transform_path))
    assert matrix.shape == (4, 4), (
        f"Expected 4x4 transform matrix, got shape {matrix.shape}"
    )

    # Sanity check: last row must be close to [0, 0, 0, 1]
    np.testing.assert_allclose(
        matrix[3, :], [0.0, 0.0, 0.0, 1.0], atol=1e-6,
        err_msg="Last row of rigid_transform.txt must be [0, 0, 0, 1]",
    )
