"""
Tier 1 analytical tests: synthetic rigid recovery.

Uses a simple 8-vertex unit cube (N=8 > NUM_FEATURES=6) so that the internal
RigidTransformer's `numElements > numFeatures` guard is satisfied.

All tests call meshmonk.rigid_register() directly — not the v0.0 stubs module.
"""

import numpy as np
import pytest

import meshmonk

# ---------------------------------------------------------------------------
# Synthetic mesh helpers
# ---------------------------------------------------------------------------

# 8-vertex unit cube — N=8, which satisfies RigidTransformer's N > 6 guard.
_CUBE_V = np.array(
    [
        [1, 1, 1],
        [1, 1, -1],
        [1, -1, 1],
        [1, -1, -1],
        [-1, 1, 1],
        [-1, 1, -1],
        [-1, -1, 1],
        [-1, -1, -1],
    ],
    dtype="float32",
)

# 12 triangles covering all 6 faces of the cube
_CUBE_F = np.array(
    [
        [0, 4, 1], [4, 5, 1],   # top  (y=+1)
        [2, 3, 6], [3, 7, 6],   # bottom (y=-1)
        [0, 2, 4], [2, 6, 4],   # front (z=+1)
        [1, 5, 3], [5, 7, 3],   # back  (z=-1)
        [0, 1, 2], [1, 3, 2],   # right (x=+1)
        [4, 6, 5], [6, 7, 5],   # left  (x=-1)
    ],
    dtype="int32",
)


def _make_rotation(axis: str, degrees: float) -> np.ndarray:
    """Return a 3x3 rotation matrix (float32) for rotation about the given axis."""
    theta = np.radians(degrees)
    c, s = np.cos(theta), np.sin(theta)
    if axis == "z":
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]], dtype="float32")
    if axis == "y":
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]], dtype="float32")
    if axis == "x":
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]], dtype="float32")
    raise ValueError(f"Unknown axis: {axis}")


def _build_transform(R: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Build 4x4 homogeneous transform from R (3,3) and t (3,)."""
    T = np.eye(4, dtype="float32")
    T[:3, :3] = R
    T[:3, 3] = t
    return T


def _apply_transform_to_features(feat: np.ndarray, R: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Apply SE(3) to (N,6) feature matrix: positions R*p+t, normals R*n."""
    out = feat.copy()
    out[:, :3] = (R @ feat[:, :3].T).T + t
    out[:, 3:] = (R @ feat[:, 3:].T).T
    return out.astype("float32")


def _rigid_register_cube(source_feat, target_feat, **kwargs):
    """Run rigid_register using raw arrays on the cube mesh."""
    flags = np.ones(source_feat.shape[0], dtype="float32")
    return meshmonk.rigid_register(
        floating_features=source_feat,
        target_features=target_feat,
        floating_faces=_CUBE_F,
        target_faces=_CUBE_F,
        floating_flags=flags,
        target_flags=flags,
        **kwargs,
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_rigid_recovery_synthetic():
    """Recover a known SE(3) from a synthetic cube.

    Applies a 15-degree rotation about Z plus a small translation to the
    floating mesh, then runs rigid_register() from identity. Asserts:
      - recovered transform matrix is close to known transform (atol=1e-3)
      - aligned_features positions close to ground truth (atol=1e-3)
      - aligned_features normals (cols 3-5) close to ground truth (atol=1e-3)
    """
    R = _make_rotation("z", 15.0)
    t = np.array([0.1, 0.2, 0.3], dtype="float32")

    feat_source = meshmonk.features_from_vertices(_CUBE_V, _CUBE_F)
    feat_target = _apply_transform_to_features(feat_source, R, t)

    result = _rigid_register_cube(feat_source, feat_target, num_iterations=80)

    T_expected = _build_transform(R, t)
    np.testing.assert_allclose(
        result.transform.matrix, T_expected, atol=1e-3,
        err_msg="Recovered transform does not match known SE(3)",
    )
    np.testing.assert_allclose(
        result.aligned_features[:, :3], feat_target[:, :3], atol=1e-3,
        err_msg="Aligned positions do not match ground truth",
    )
    np.testing.assert_allclose(
        result.aligned_features[:, 3:], feat_target[:, 3:], atol=1e-3,
        err_msg="Aligned normals do not match ground truth",
    )


def test_self_consistency():
    """Register A→B then B→A; assert T1.compose(T2) ~ Identity (atol=1e-3)."""
    R = _make_rotation("y", 10.0)
    t = np.array([0.05, 0.1, 0.15], dtype="float32")

    feat_a = meshmonk.features_from_vertices(_CUBE_V, _CUBE_F)
    feat_b = _apply_transform_to_features(feat_a, R, t)

    res_ab = _rigid_register_cube(feat_a, feat_b, num_iterations=80)
    res_ba = _rigid_register_cube(feat_b, feat_a, num_iterations=80)

    composed = res_ab.transform.compose(res_ba.transform)
    np.testing.assert_allclose(
        composed.matrix, np.eye(4, dtype="float32"), atol=1e-3,
        err_msg="T(A→B).compose(T(B→A)) is not close to Identity",
    )


def test_round_trip():
    """Apply known transform, register, compose result with inverse → Identity."""
    R = _make_rotation("x", 20.0)
    t = np.array([0.2, -0.1, 0.1], dtype="float32")

    feat_a = meshmonk.features_from_vertices(_CUBE_V, _CUBE_F)
    feat_b = _apply_transform_to_features(feat_a, R, t)

    result = _rigid_register_cube(feat_a, feat_b, num_iterations=80)

    T_known = _build_transform(R, t)
    T_known_rt = meshmonk.RigidTransform()
    # Build inverse by direct matrix math (T_known^{-1})
    # For SE(3): [R|t]^{-1} = [R^T | -R^T*t]
    R_inv = R.T
    t_inv = -(R_inv @ t)
    T_inv = _build_transform(R_inv, t_inv)
    # Create RigidTransform from inverse using compute_rigid_transform trick:
    # compose recovered transform with the known inverse
    # We build an identity RigidTransform and "wrap" the matrix manually
    # Since RigidTransform.matrix has no setter, we use matrix math:
    T_composed = result.transform.matrix @ T_inv
    np.testing.assert_allclose(
        T_composed, np.eye(4, dtype="float32"), atol=1e-3,
        err_msg="Recovered transform composed with known inverse is not close to Identity",
    )


def test_correspondence_sanity():
    """Identical meshes → correspondences map each vertex to itself (atol=1e-4)."""
    feat = meshmonk.features_from_vertices(_CUBE_V, _CUBE_F)
    flags = np.ones(feat.shape[0], dtype="float32")

    corr_feat, corr_flags = meshmonk.compute_correspondences(feat, feat, flags, flags)

    np.testing.assert_allclose(
        corr_feat, feat, atol=1e-4,
        err_msg="Correspondences on identical meshes should be identity",
    )


def test_inlier_weight_sanity():
    """Identical meshes → inlier weights near 1.0 (atol=1e-4)."""
    feat = meshmonk.features_from_vertices(_CUBE_V, _CUBE_F)
    flags = np.ones(feat.shape[0], dtype="float32")

    corr_feat, corr_flags = meshmonk.compute_correspondences(feat, feat, flags, flags)
    weights = meshmonk.compute_inlier_weights(feat, corr_feat, corr_flags)

    np.testing.assert_allclose(
        weights, np.ones(feat.shape[0], dtype="float32"), atol=1e-4,
        err_msg="Inlier weights on identical meshes should be near 1.0",
    )
