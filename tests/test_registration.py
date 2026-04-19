"""Tier 1/2 analytical tests for nonrigid and pyramid registration."""

import numpy as np
import meshmonk


def _make_test_features():
    """Create a simple mesh for testing — a cube with computed normals.

    Returns (features, faces, flags) ready for Pattern B registration.
    """
    # 8 vertices of a unit cube
    vertices = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ],
        dtype=np.float32,
    )

    faces = np.array(
        [
            [0, 1, 2],
            [0, 2, 3],  # bottom
            [4, 6, 5],
            [4, 7, 6],  # top
            [0, 4, 5],
            [0, 5, 1],  # front
            [2, 6, 7],
            [2, 7, 3],  # back
            [0, 3, 7],
            [0, 7, 4],  # left
            [1, 5, 6],
            [1, 6, 2],  # right
        ],
        dtype=np.int32,
    )

    flags = np.ones(8, dtype=np.float32)

    # Compute proper normals via meshmonk
    normals = np.asarray(meshmonk.compute_normals(vertices, faces), dtype=np.float32)
    features = np.empty((8, 6), dtype=np.float32)
    features[:, :3] = vertices
    features[:, 3:] = normals

    return features, faces, flags


class TestNonrigidSelfRegistration:
    """Tier 1: Nonrigid registration of a mesh against itself should produce near-zero displacement."""

    def test_self_registration_near_identity(self):
        """Registering a mesh to itself should produce minimal displacement."""
        features, faces, flags = _make_test_features()

        result = meshmonk.nonrigid_register(
            floating_features=features.copy(),
            target_features=features.copy(),
            floating_faces=faces,
            target_faces=faces,
            floating_flags=flags,
            target_flags=flags,
            num_iterations=10,
        )

        # Displacement should be very small for self-registration
        disp = result.displacement_field
        assert disp.shape == (8, 3)
        rmse = np.sqrt(np.mean(disp**2))
        assert rmse < 1.0, f"Self-registration displacement RMSE {rmse:.4f} exceeds 1.0"

    def test_result_shapes(self):
        """Verify nonrigid result has all expected fields with correct shapes."""
        features, faces, flags = _make_test_features()

        result = meshmonk.nonrigid_register(
            floating_features=features.copy(),
            target_features=features.copy(),
            floating_faces=faces,
            target_faces=faces,
            floating_flags=flags,
            target_flags=flags,
            num_iterations=10,
        )

        assert result.aligned_features.shape == (8, 6)
        assert result.displacement_field.shape == (8, 3)
        assert result.final_inlier_weights.shape == (8,)
        assert result.iterations_run == 10


class TestPyramidSelfRegistration:
    """Tier 1: Pyramid registration of a mesh against itself."""

    def test_self_registration_near_identity(self):
        """Registering a mesh to itself should produce minimal displacement."""
        features, faces, flags = _make_test_features()

        result = meshmonk.pyramid_register(
            floating_features=features.copy(),
            target_features=features.copy(),
            floating_faces=faces,
            target_faces=faces,
            floating_flags=flags,
            target_flags=flags,
            num_iterations=6,
            num_pyramid_layers=1,
        )

        disp = result.displacement_field
        assert disp.shape == (8, 3)
        rmse = np.sqrt(np.mean(disp**2))
        assert rmse < 1.0, f"Self-registration displacement RMSE {rmse:.4f} exceeds 1.0"

    def test_result_shapes_and_layers(self):
        """Verify pyramid result has all expected fields."""
        features, faces, flags = _make_test_features()

        result = meshmonk.pyramid_register(
            floating_features=features.copy(),
            target_features=features.copy(),
            floating_faces=faces,
            target_faces=faces,
            floating_flags=flags,
            target_flags=flags,
            num_iterations=6,
            num_pyramid_layers=2,
        )

        assert result.aligned_features.shape == (8, 6)
        assert result.displacement_field.shape == (8, 3)
        assert result.final_inlier_weights.shape == (8,)
        assert len(result.per_layer_iterations) == 2
        assert sum(result.per_layer_iterations) > 0

    def test_single_layer_matches_nonrigid(self):
        """Pyramid with 1 layer should produce similar results to plain nonrigid."""
        features, faces, flags = _make_test_features()

        nr_result = meshmonk.nonrigid_register(
            floating_features=features.copy(),
            target_features=features.copy(),
            floating_faces=faces,
            target_faces=faces,
            floating_flags=flags,
            target_flags=flags,
            num_iterations=6,
        )

        pyr_result = meshmonk.pyramid_register(
            floating_features=features.copy(),
            target_features=features.copy(),
            floating_faces=faces,
            target_faces=faces,
            floating_flags=flags,
            target_flags=flags,
            num_iterations=6,
            num_pyramid_layers=1,
        )

        # Both should produce small displacements for self-registration
        nr_rmse = np.sqrt(np.mean(nr_result.displacement_field**2))
        pyr_rmse = np.sqrt(np.mean(pyr_result.displacement_field**2))
        # Both should be small (self-registration)
        assert nr_rmse < 1.0
        assert pyr_rmse < 1.0
        assert abs(nr_rmse - pyr_rmse) < 0.5


class TestNonrigidWithTranslation:
    """Tier 2: Nonrigid should reduce displacement when given a translated mesh."""

    def test_translated_mesh_converges(self):
        """A small translation should be (partially) recovered by nonrigid registration."""
        features, faces, flags = _make_test_features()

        # Create a translated copy
        translated = features.copy()
        translated[:, 0] += 0.3  # shift X by 0.3 units

        result = meshmonk.nonrigid_register(
            floating_features=translated,
            target_features=features.copy(),
            floating_faces=faces,
            target_faces=faces,
            floating_flags=flags,
            target_flags=flags,
            num_iterations=20,
        )

        # The aligned vertices should be closer to the target than the original
        original_dist = np.sqrt(np.mean((translated[:, :3] - features[:, :3]) ** 2))
        aligned_dist = np.sqrt(
            np.mean((result.aligned_vertices - features[:, :3]) ** 2)
        )
        assert aligned_dist < original_dist, (
            f"Nonrigid registration should reduce distance to target. "
            f"Original RMSE: {original_dist:.4f}, Aligned RMSE: {aligned_dist:.4f}"
        )
        assert aligned_dist < 0.5 * original_dist, (
            f"Nonrigid registration should substantially reduce distance to target "
            f"(expected < 50% of original). "
            f"Original RMSE: {original_dist:.4f}, Aligned RMSE: {aligned_dist:.4f}"
        )
