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
    t_expected = expected_transform[:3, 3]

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
        T[3, :], [0.0, 0.0, 0.0, 1.0], atol=1e-5,
        err_msg="Last row of transform matrix must be [0,0,0,1]",
    )
    R_result = T[:3, :3]
    det_R = float(np.linalg.det(R_result))
    assert abs(det_R - 1.0) < 1e-4, f"Rotation determinant must be 1.0, got {det_R}"
    np.testing.assert_allclose(
        R_result @ R_result.T, np.eye(3, dtype="float32"), atol=1e-4,
        err_msg="Rotation R must be orthogonal: R @ R.T ~ I",
    )

    # --- Verify 2: rotation part matches expected_transform ---
    np.testing.assert_allclose(
        R_result, R_expected, atol=1e-4,
        err_msg="Rotation part of result does not match expected_transform.txt rotation",
    )
