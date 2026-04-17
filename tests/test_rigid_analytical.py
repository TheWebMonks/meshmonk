"""
Tier 1 analytical test: synthetic rigid recovery.

This test is INTENTIONALLY RED at v0.0 — the stub raises NotImplementedError.
It ships as a skeleton; v0.1 implements recover_rigid and makes it green.

v0.0 contract: pytest exits with code 1 (exactly this test errors).
"""

import numpy as np
import pytest
from tests.stubs.rigid import recover_rigid


def _apply_rigid(vertices: np.ndarray, R: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Apply rotation R (3,3) and translation t (3,) to (N,3) vertices."""
    return (R @ vertices.T).T + t


@pytest.mark.xfail(raises=NotImplementedError, strict=True)
def test_synthetic_rigid_recovery():
    """Recover a known SE(3) from a synthetic 20-vertex point cloud.

    Build a source, apply a known rotation + translation, then call
    recover_rigid(src, tgt) and assert the recovered transform matches.

    EXPECTED TO FAIL at v0.0 with NotImplementedError.
    First green in v0.1.
    """
    rng = np.random.default_rng(42)
    src = rng.standard_normal((20, 3)).astype(np.float64)

    # Known SE(3): 30-degree rotation around Z + translation (1, 2, 3)
    theta = np.radians(30.0)
    R_known = np.array(
        [
            [np.cos(theta), -np.sin(theta), 0.0],
            [np.sin(theta), np.cos(theta), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    t_known = np.array([1.0, 2.0, 3.0])
    tgt = _apply_rigid(src, R_known, t_known)

    # recover_rigid raises NotImplementedError at v0.0
    T_recovered = recover_rigid(src, tgt)  # (4, 4) homogeneous transform

    R_recovered = T_recovered[:3, :3]
    t_recovered = T_recovered[:3, 3]

    np.testing.assert_allclose(R_recovered, R_known, atol=1e-4)
    np.testing.assert_allclose(t_recovered, t_known, atol=1e-4)
