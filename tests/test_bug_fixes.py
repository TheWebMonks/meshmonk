"""
Tests for critical bug fixes.

Bug 1: transform.hpp inverse() — squaredNorm vs norm regression
Bug 2: pyramid_registration() num_pyramid_layers=0 causes division-by-zero UB
Bug 3: pyramid_registration() iters_per_layer < 2 causes internal division-by-zero

Notes on test strategy:
  - RigidTransform.matrix is read-only from Python (def_prop_ro), so we cannot
    assign a matrix directly. Instead we obtain RigidTransform objects via
    compute_rigid_transform(use_scaling=True) with synthetic point clouds.
  - pyramid_registration() in Python raises MeshMonkError (not returns Expected),
    so we use pytest.raises().
"""

from __future__ import annotations

import numpy as np
import pytest

import meshmonk
from meshmonk import MeshMonkError


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# 8-vertex unit cube
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

_CUBE_F = np.array(
    [
        [0, 4, 1], [4, 5, 1],
        [2, 3, 6], [3, 7, 6],
        [0, 2, 4], [2, 6, 4],
        [1, 5, 3], [5, 7, 3],
        [0, 1, 2], [1, 3, 2],
        [4, 6, 5], [6, 7, 5],
    ],
    dtype="int32",
)


def _make_cube_features():
    """Return (N, 6) feature matrix for the unit cube."""
    return meshmonk.features_from_vertices(_CUBE_V, _CUBE_F)


def _make_cube_mesh():
    """Return (features, faces, flags) for the unit cube."""
    features = _make_cube_features()
    flags = np.ones(_CUBE_V.shape[0], dtype="float32")
    return features, _CUBE_F, flags


def _get_rigid_transform_s1():
    """Return a RigidTransform with known pure rotation (s=1) via compute_rigid_transform.

    Applies 30-degree rotation about Z to the cube, then uses compute_rigid_transform
    with use_scaling=False to recover a unit-scale transform.
    """
    theta = np.radians(30)
    R = np.array(
        [
            [np.cos(theta), -np.sin(theta), 0],
            [np.sin(theta),  np.cos(theta), 0],
            [0,              0,             1],
        ],
        dtype="float32",
    )
    t = np.array([1.0, 2.0, 3.0], dtype="float32")
    features = _make_cube_features()

    # Build correspondence: apply known transform to floating features
    corr = features.copy()
    corr[:, :3] = (R @ features[:, :3].T).T + t
    corr[:, 3:] = (R @ features[:, 3:].T).T

    weights = np.ones(8, dtype="float32")
    T = meshmonk.compute_rigid_transform(features, corr, weights, use_scaling=False)
    return T, R, t


def _get_rigid_transform_scaled(scale: float):
    """Return a RigidTransform with known scale factor via compute_rigid_transform(use_scaling=True).

    Applies (scale * R + t) to the cube. Uses a well-conditioned setup so that
    compute_rigid_transform with use_scaling=True returns the correct sR matrix.
    """
    theta = np.radians(30)
    R_pure = np.array(
        [
            [np.cos(theta), -np.sin(theta), 0],
            [np.sin(theta),  np.cos(theta), 0],
            [0,              0,             1],
        ],
        dtype="float32",
    )
    sR = scale * R_pure
    t = np.array([5.0, -3.0, 2.0], dtype="float32")
    features = _make_cube_features()

    # Build perfect correspondences under the similarity transform
    corr = features.copy()
    corr[:, :3] = (sR @ features[:, :3].T).T + t
    corr[:, 3:] = (R_pure @ features[:, 3:].T).T  # normals use pure rotation

    weights = np.ones(8, dtype="float32")
    T = meshmonk.compute_rigid_transform(features, corr, weights, use_scaling=True)
    return T, sR, t


# ---------------------------------------------------------------------------
# Bug 1: RigidTransform.inverse() squaredNorm fix
# ---------------------------------------------------------------------------


class TestTransformInverse:
    """Bug 1: inverse() must use squaredNorm (divide by s^2) not norm (divide by s)."""

    def test_pure_rotation_inverse_is_identity(self):
        """T.matrix @ T.inverse().matrix == Identity for a pure rotation (s=1).

        With the norm() bug, this would also fail due to dividing by s=1 (no
        observable difference), but a correct identity ensures the math is right.
        """
        T, R, t = _get_rigid_transform_s1()

        T_inv = T.inverse()
        product = T.matrix @ T_inv.matrix
        np.testing.assert_allclose(
            product, np.eye(4, dtype="float32"), atol=1e-4,
            err_msg="T * T^-1 must be Identity for pure rotation (s=1)",
        )

    def test_similarity_inverse_is_identity(self):
        """T.matrix @ T.inverse().matrix == Identity for a similarity transform (s=2).

        This is the critical test that catches the norm() vs squaredNorm() bug:
        - With norm():        T_inv upper-left = Rt / s  = R_pure^T  (missing 1/s)
        - With squaredNorm(): T_inv upper-left = Rt / s2 = R_pure^T / s (correct)
        So T * T_inv has scale factor s in the upper-left block instead of I.
        """
        scale = 2.0
        T, sR, t = _get_rigid_transform_scaled(scale)

        # Verify the transform has the right scale embedded
        s_actual = float(np.linalg.norm(T.matrix[:3, 0]))  # col-0 norm = s
        np.testing.assert_allclose(s_actual, scale, atol=0.05,
            err_msg="compute_rigid_transform did not return expected scale")

        T_inv = T.inverse()
        product = T.matrix @ T_inv.matrix
        np.testing.assert_allclose(
            product, np.eye(4, dtype="float32"), atol=1e-4,
            err_msg=(
                "T * T^-1 must be Identity for similarity transform (scale=2). "
                "Failure here means inverse() divides by s instead of s^2."
            ),
        )

    def test_inverse_apply_roundtrip_scaled(self):
        """Applying T then T.inverse() to positions recovers original positions.

        Uses scale=3 to maximally stress the squaredNorm fix.
        This is an end-to-end test: build transform -> apply -> invert -> recover.
        """
        scale = 3.0
        T, sR, t = _get_rigid_transform_scaled(scale)
        features = _make_cube_features()

        # Apply T then T.inverse()
        transformed = T.apply(features)
        recovered = T.inverse().apply(transformed)

        np.testing.assert_allclose(
            recovered[:, :3], features[:, :3], atol=1e-3,
            err_msg=(
                "Applying T then T^-1 must recover original positions. "
                "Failure here means inverse() uses wrong scale factor."
            ),
        )


# ---------------------------------------------------------------------------
# Bug 2 + 3: Pyramid registration validation
# ---------------------------------------------------------------------------


class TestPyramidValidation:
    """Bugs 2 & 3: pyramid_registration() must validate pyramid params before running."""

    def _call_pyramid(self, num_iterations: int, num_pyramid_layers: int):
        """Call _pyramid_registration directly with minimal valid mesh data."""
        from meshmonk._meshmonk_core import (  # noqa: PLC0415
            pyramid_registration as _pyramid_registration,
            PyramidParams,
        )

        features, faces, flags = _make_cube_mesh()
        params = PyramidParams()
        params.num_iterations = num_iterations
        params.num_pyramid_layers = num_pyramid_layers

        return _pyramid_registration(
            features, features,
            faces, faces,
            flags, flags,
            params,
        )

    def test_zero_pyramid_layers_raises_degenerate_error(self):
        """num_pyramid_layers=0 must raise MeshMonkError with DegenerateInput message.

        Bug 2: without the fix, divides by zero at iters_per_layer computation.
        Currently the C++ logs a warning but proceeds — no exception is raised.
        """
        with pytest.raises(MeshMonkError, match="DegenerateInput"):
            self._call_pyramid(num_iterations=100, num_pyramid_layers=0)

    def test_too_few_iters_per_layer_raises_degenerate_error(self):
        """iters_per_layer < 2 must raise MeshMonkError with DegenerateInput message.

        Bug 3: 2 iterations / 3 layers = round(2/3) = 1 iteration per layer,
        which causes division-by-zero inside NonrigidRegistration:
          exp(log(end/start) / (iters - 1)) with iters=1.
        """
        with pytest.raises(MeshMonkError, match="DegenerateInput"):
            self._call_pyramid(num_iterations=2, num_pyramid_layers=3)

    def test_one_iteration_one_layer_raises_degenerate_error(self):
        """1 iteration / 1 layer = 1 per layer — must be rejected (need >= 2)."""
        with pytest.raises(MeshMonkError, match="DegenerateInput"):
            self._call_pyramid(num_iterations=1, num_pyramid_layers=1)

    def test_valid_params_do_not_raise_degenerate_error(self):
        """2 iterations / 1 layer = 2 per layer — must NOT raise DegenerateInput.

        This guards against over-eager validation.
        Any other error (e.g. InsufficientInliers) is acceptable.
        """
        try:
            self._call_pyramid(num_iterations=2, num_pyramid_layers=1)
        except MeshMonkError as e:
            assert "DegenerateInput" not in str(e), (
                "pyramid_registration with 2 iterations / 1 layer (2 per layer) "
                "must not raise DegenerateInput from the new validation"
            )


# ---------------------------------------------------------------------------
# Bug 4: PyramidNonrigidRegistration downsample target guard copy-paste bug
# ---------------------------------------------------------------------------


class TestPyramidDownsampleTargetGuard:
    """Bug 4: line 54 of PyramidNonrigidRegistration.cpp had a copy-paste bug.

    The guard for _downsampleTargetStart checked _downsampleFloatStart instead:

        WRONG (bug):
            if ((_downsampleTargetStart < 0.0f) || (_downsampleFloatStart >= 100.0f))

        CORRECT (fix):
            if ((_downsampleTargetStart < 0.0f) || (_downsampleTargetStart >= 100.0f))

    This means an out-of-range _downsampleTargetStart (>= 100) would not be
    clamped to 90.0f when _downsampleFloatStart is valid (< 100).

    The clamping is internal state — not directly observable from Python (C++
    stdout is not capturable through nanobind with gil_scoped_release, and the
    downsampler silently falls back to keeping all vertices with an invalid ratio).

    We therefore verify the fix via source-code inspection, which directly
    captures the exact bug and its correction.
    """

    _SOURCE_FILE = "/workspace/library/src/PyramidNonrigidRegistration.cpp"

    def test_target_guard_checks_target_variable(self):
        """Line 54: both conditions in the _downsampleTargetStart guard must
        reference _downsampleTargetStart, not _downsampleFloatStart.

        Fails before the fix (bug present) and passes after the fix.
        """
        with open(self._SOURCE_FILE) as f:
            source = f.read()

        # The fixed line: both conditions must check _downsampleTargetStart
        correct_guard = (
            "if ((_downsampleTargetStart < 0.0f) || (_downsampleTargetStart >= 100.0f))"
        )
        # The buggy line: second condition incorrectly checks _downsampleFloatStart
        buggy_guard = (
            "if ((_downsampleTargetStart < 0.0f) || (_downsampleFloatStart >= 100.0f))"
        )

        assert buggy_guard not in source, (
            "Copy-paste bug found in PyramidNonrigidRegistration.cpp line 54: "
            "the _downsampleTargetStart guard checks _downsampleFloatStart instead "
            "of _downsampleTargetStart in its second condition. "
            f"Remove: {buggy_guard!r}"
        )
        assert correct_guard in source, (
            "Expected the corrected _downsampleTargetStart guard to be present. "
            f"Expected: {correct_guard!r}"
        )


# ---------------------------------------------------------------------------
# Target bounding box degeneracy validation
# ---------------------------------------------------------------------------


class TestTargetBboxValidation:
    """Verify target mesh with degenerate bounding box raises MeshMonkError."""

    def _make_degenerate_target(self):
        """Target where all vertices are at the same point."""
        n = 10
        # Normal floating mesh
        float_verts = np.random.randn(n, 3).astype(np.float32) * 10
        float_normals = np.random.randn(n, 3).astype(np.float32)
        float_feat = np.hstack([float_verts, float_normals])
        float_faces = np.array([[0, 1, 2]], dtype=np.int32)
        float_flags = np.ones(n, dtype=np.float32)

        # Degenerate target: all vertices at same point
        target_verts = np.ones((n, 3), dtype=np.float32) * 5.0
        target_normals = np.random.randn(n, 3).astype(np.float32)
        target_feat = np.hstack([target_verts, target_normals])
        target_faces = np.array([[0, 1, 2]], dtype=np.int32)
        target_flags = np.ones(n, dtype=np.float32)

        return (float_feat, target_feat, float_faces, target_faces,
                float_flags, target_flags)

    def test_rigid_degenerate_target_raises(self):
        args = self._make_degenerate_target()
        with pytest.raises(meshmonk.MeshMonkError):
            meshmonk.rigid_register(
                floating_features=args[0], target_features=args[1],
                floating_faces=args[2], target_faces=args[3],
                floating_flags=args[4], target_flags=args[5],
            )

    def test_nonrigid_degenerate_target_raises(self):
        args = self._make_degenerate_target()
        with pytest.raises(meshmonk.MeshMonkError):
            meshmonk.nonrigid_register(
                floating_features=args[0], target_features=args[1],
                floating_faces=args[2], target_faces=args[3],
                floating_flags=args[4], target_flags=args[5],
            )

    def test_pyramid_degenerate_target_raises(self):
        args = self._make_degenerate_target()
        with pytest.raises(meshmonk.MeshMonkError):
            meshmonk.pyramid_register(
                floating_features=args[0], target_features=args[1],
                floating_faces=args[2], target_faces=args[3],
                floating_flags=args[4], target_flags=args[5],
            )
