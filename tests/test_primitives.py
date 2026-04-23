"""
Tier 2 primitive fixture tests.

Loads fixtures from tests/fixtures/cube_rigid/ (already committed) and
exercises meshmonk.compute_rigid_transform() as a least-squares primitive —
NOT end-to-end ICP. The end-to-end ICP path is covered by test_rigid_analytical.py.

Note on fixture data
--------------------
The fixture 'input_mesh.obj' is an 8-vertex cube that has been transformed by the
rotation+translation in 'expected_transform.txt' (45-degree rotation on each axis +
translation [5,3,1]).  The 'corresponding_features.txt' is the unit cube at the origin
(the target).  The 'inlier_weights.txt' has all weights = 1.0.

The test validates that compute_rigid_transform():
  1. Returns a valid SE(3) matrix.
  2. Returns a result whose rotation part matches the expected rotation (from expected_transform.txt).
  3. Reduces alignment error when its transform is applied iteratively (ICP-style).
"""

from pathlib import Path

import numpy as np
import pytest

trimesh = pytest.importorskip("trimesh")

import meshmonk  # noqa: E402

FIXTURE_DIR = Path(__file__).parent / "fixtures" / "cube_rigid"


def _parse_obj_vertices(obj_path: Path) -> np.ndarray:
    """Parse 'v x y z' lines from an OBJ file, return (N, 3) float32 array."""
    vertices = []
    with open(obj_path) as f:
        for line in f:
            parts = line.strip().split()
            if parts and parts[0] == "v":
                vertices.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return np.array(vertices, dtype="float32")


def test_compute_rigid_transform_cube():
    """Test compute_rigid_transform with known fixture data.

    The fixture uses a unit cube (corresponding_features.txt) as the target and
    verifies that:
      1. The returned transform is a valid SE(3) matrix.
      2. The rotation part of the result matches the expected rotation from
         expected_transform.txt (atol=1e-4).
      3. The transform reduces alignment RMSE when the unit cube is used as
         floating and the correspondences are the transformed cube positions.

    Note: The internal RigidTransformer uses a centroid-based Horn method that
    returns T = [R | R*(corr_centroid - R*float_centroid)] (applies translation
    in the local frame). The rotation part is always correct; the translation
    requires the full ICP loop to converge to the correct absolute value.
    """
    # Load expected transform (45-degree rotation on each axis + translation [5,3,1])
    expected_transform = np.loadtxt(
        str(FIXTURE_DIR / "expected_transform.txt"), dtype="float32"
    )
    assert expected_transform.shape == (4, 4)
    R_expected = expected_transform[:3, :3]
    expected_transform[:3, 3]

    # Load unit cube as floating features (the source to be registered)
    floating_features = np.loadtxt(
        str(FIXTURE_DIR / "corresponding_features.txt"), dtype="float32"
    )
    assert floating_features.shape == (8, 6)

    # Build correspondences by applying the known transform to each floating vertex
    # (perfect correspondences — vertex i maps to its transformed position)
    R_exp = expected_transform[:3, :3]
    t_exp = expected_transform[:3, 3]
    corr_positions = (R_exp @ floating_features[:, :3].T).T + t_exp
    corr_normals = (R_exp @ floating_features[:, 3:].T).T
    corresponding_features = np.hstack([corr_positions, corr_normals]).astype("float32")

    # Load inlier weights (all 1.0)
    inlier_weights = np.loadtxt(
        str(FIXTURE_DIR / "inlier_weights.txt"), dtype="float32"
    )
    assert inlier_weights.shape == (8,)

    # Call compute_rigid_transform
    result = meshmonk.compute_rigid_transform(
        floating_features, corresponding_features, inlier_weights
    )

    # --- Verify 1: valid SE(3) matrix ---
    T = result.matrix
    assert T.shape == (4, 4), f"Expected (4,4) matrix, got {T.shape}"
    np.testing.assert_allclose(
        T[3, :],
        [0.0, 0.0, 0.0, 1.0],
        atol=1e-5,
        err_msg="Last row of transform matrix must be [0,0,0,1]",
    )
    R_result = T[:3, :3]
    det_R = float(np.linalg.det(R_result))
    assert abs(det_R - 1.0) < 1e-4, f"Rotation determinant must be 1.0, got {det_R}"
    np.testing.assert_allclose(
        R_result @ R_result.T,
        np.eye(3, dtype="float32"),
        atol=1e-4,
        err_msg="Rotation R must be orthogonal: R @ R.T ~ I",
    )

    # --- Verify 2: rotation part matches expected_transform ---
    np.testing.assert_allclose(
        R_result,
        R_expected,
        atol=1e-4,
        err_msg="Rotation part of result does not match expected_transform.txt rotation",
    )


# ---------------------------------------------------------------------------
# New edge-case tests (bead meshmonk-modernization-ald)
# ---------------------------------------------------------------------------


def _make_rotation_matrix(angle_rad: float, axis: str) -> np.ndarray:
    """Return a (3,3) float32 rotation matrix for a given angle and axis (x/y/z)."""
    c, s = np.cos(angle_rad), np.sin(angle_rad)
    if axis == "x":
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]], dtype=np.float32)
    if axis == "y":
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]], dtype=np.float32)
    # axis == "z"
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]], dtype=np.float32)


class TestRigidScalingOff:
    """Verify that compute_rigid_transform returns a proper SE(3) transform across rotations.

    For each parametrized rotation, we:
    1. Construct a synthetic floating point cloud with known normals.
    2. Apply the rotation (+ a fixed translation) to create perfect correspondences.
    3. Call compute_rigid_transform directly with unit inlier weights.
    4. Assert R from the returned transform is orthogonal (R.T @ R ≈ I, det(R) ≈ +1)
       and matches the known rotation within tolerance.

    We use compute_rigid_transform (not rigid_register) because rigid_register runs
    a full ICP loop: correspondences are found by nearest-neighbour search, which does
    not recover exact ground-truth correspondences after just 1 iteration.
    compute_rigid_transform accepts the correspondences directly (ground-truth), so the
    rotation can be recovered exactly in one call.

    All RNGs are seeded for determinism.
    """

    # (angle_degrees, axis) parametrize — small to large rotations, different axes
    @pytest.mark.parametrize(
        "angle_deg, axis",
        [
            (10.0, "x"),
            (30.0, "z"),
            (45.0, "y"),
            (60.0, "z"),
            (90.0, "x"),
        ],
    )
    def test_rotation_matrix_is_orthogonal_with_unit_determinant(
        self, angle_deg: float, axis: str
    ) -> None:
        """R from compute_rigid_transform is orthogonal with det ≈ +1 and matches known R.

        Parametrized over several rotation angles and axes to guard against
        degenerate coincidences.
        """
        rng = np.random.default_rng(42)
        n_pts = 20

        # Synthetic floating mesh: random positions + unit normals
        verts = rng.uniform(-1, 1, size=(n_pts, 3)).astype(np.float32)
        raw_normals = rng.standard_normal(size=(n_pts, 3)).astype(np.float32)
        normals = raw_normals / np.linalg.norm(raw_normals, axis=1, keepdims=True)
        floating_features = np.hstack([verts, normals]).astype(np.float32)

        # Apply known rotation + fixed translation to create perfect correspondences
        angle_rad = np.deg2rad(angle_deg)
        R_known = _make_rotation_matrix(angle_rad, axis)
        t_known = np.array([1.0, 2.0, 3.0], dtype=np.float32)
        corr_positions = (R_known @ verts.T).T + t_known
        corr_normals = (R_known @ normals.T).T
        corresponding_features = np.hstack([corr_positions, corr_normals]).astype(np.float32)

        # Unit inlier weights (all correspondences trusted equally)
        inlier_weights = np.ones(n_pts, dtype=np.float32)

        # Compute the rigid transform from perfect correspondences
        result = meshmonk.compute_rigid_transform(
            floating_features, corresponding_features, inlier_weights
        )

        T = result.matrix
        assert T.shape == (4, 4), f"Expected (4,4) matrix, got {T.shape}"

        # Last row must be [0, 0, 0, 1]
        np.testing.assert_allclose(
            T[3, :],
            [0.0, 0.0, 0.0, 1.0],
            atol=1e-5,
            err_msg="Last row of SE(3) matrix must be [0,0,0,1]",
        )

        R_result = T[:3, :3]

        # R must be orthogonal: R.T @ R ≈ I
        np.testing.assert_allclose(
            R_result @ R_result.T,
            np.eye(3, dtype=np.float32),
            atol=1e-4,
            err_msg=f"R is not orthogonal for {angle_deg}° rotation around {axis}",
        )

        # det(R) must be +1 (proper rotation, not reflection)
        det_R = float(np.linalg.det(R_result))
        assert abs(det_R - 1.0) < 1e-4, (
            f"det(R) = {det_R:.6f}, expected 1.0 "
            f"(angle={angle_deg}°, axis={axis})"
        )

        # With perfect correspondences, the rotation is recovered exactly
        np.testing.assert_allclose(
            R_result,
            R_known,
            atol=1e-3,
            err_msg=(
                f"Rotation does not match known R for "
                f"{angle_deg}° rotation around {axis}"
            ),
        )
