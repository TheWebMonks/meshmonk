"""
Integration tests for the Tier 3.5 legacy baseline golden fixtures.

These tests verify that the ingested .npz goldens from demo/rigid_output.obj
and demo/pyramid_output.obj satisfy the pinned schema and are genuinely distinct
meshes (mean per-vertex displacement > 0.1 mm).
"""

from pathlib import Path

import numpy as np
import pytest

GOLDEN_DIR = Path(__file__).parent / "legacy_baseline"


@pytest.fixture
def rigid_npz():
    path = GOLDEN_DIR / "rigid_output.npz"
    if not path.exists():
        pytest.skip(f"Golden file not yet generated: {path}")
    return np.load(path)


@pytest.fixture
def pyramid_npz():
    path = GOLDEN_DIR / "pyramid_output.npz"
    if not path.exists():
        pytest.skip(f"Golden file not yet generated: {path}")
    return np.load(path)


class TestLegacyBaselineSchema:
    """Verify pinned schema: vertices (N,3) float64, faces (M,3) int32."""

    def test_rigid_npz_exists(self):
        assert (
            GOLDEN_DIR / "rigid_output.npz"
        ).exists(), "rigid_output.npz not found — run conversion step"

    def test_pyramid_npz_exists(self):
        assert (
            GOLDEN_DIR / "pyramid_output.npz"
        ).exists(), "pyramid_output.npz not found — run conversion step"

    def test_rigid_vertices_shape(self, rigid_npz):
        assert "vertices" in rigid_npz, "rigid_output.npz missing 'vertices' key"
        assert rigid_npz["vertices"].ndim == 2
        assert (
            rigid_npz["vertices"].shape[1] == 3
        ), f"Expected (N,3) vertices, got {rigid_npz['vertices'].shape}"

    def test_rigid_vertices_dtype(self, rigid_npz):
        assert (
            rigid_npz["vertices"].dtype == np.float64
        ), f"Expected float64 vertices, got {rigid_npz['vertices'].dtype}"

    def test_rigid_faces_shape(self, rigid_npz):
        assert "faces" in rigid_npz, "rigid_output.npz missing 'faces' key"
        assert rigid_npz["faces"].ndim == 2
        assert (
            rigid_npz["faces"].shape[1] == 3
        ), f"Expected (M,3) faces, got {rigid_npz['faces'].shape}"

    def test_rigid_faces_dtype(self, rigid_npz):
        assert (
            rigid_npz["faces"].dtype == np.int32
        ), f"Expected int32 faces, got {rigid_npz['faces'].dtype}"

    def test_pyramid_vertices_shape(self, pyramid_npz):
        assert "vertices" in pyramid_npz
        assert pyramid_npz["vertices"].ndim == 2
        assert pyramid_npz["vertices"].shape[1] == 3

    def test_pyramid_vertices_dtype(self, pyramid_npz):
        assert pyramid_npz["vertices"].dtype == np.float64

    def test_pyramid_faces_shape(self, pyramid_npz):
        assert "faces" in pyramid_npz
        assert pyramid_npz["faces"].ndim == 2
        assert pyramid_npz["faces"].shape[1] == 3

    def test_pyramid_faces_dtype(self, pyramid_npz):
        assert pyramid_npz["faces"].dtype == np.int32


class TestLegacyBaselineDistinct:
    """Verify rigid and pyramid meshes are genuinely different."""

    def test_same_vertex_count(self, rigid_npz, pyramid_npz):
        """Both outputs should map the same source mesh — vertex count must match."""
        n_rigid = rigid_npz["vertices"].shape[0]
        n_pyramid = pyramid_npz["vertices"].shape[0]
        assert (
            n_rigid == n_pyramid
        ), f"Vertex count mismatch: rigid={n_rigid}, pyramid={n_pyramid}"

    def test_mean_displacement_exceeds_threshold(self, rigid_npz, pyramid_npz):
        """Mean per-vertex displacement between rigid and pyramid must exceed 0.1 mm.

        Guards against the two goldens being accidental copies of each other.
        """
        rigid_verts = rigid_npz["vertices"]
        pyramid_verts = pyramid_npz["vertices"]
        displacements = np.linalg.norm(rigid_verts - pyramid_verts, axis=1)
        mean_displacement = float(np.mean(displacements))
        assert mean_displacement > 0.1, (
            f"Mean per-vertex displacement {mean_displacement:.4f} mm <= 0.1 mm; "
            "rigid and pyramid goldens may be identical"
        )


class TestLegacyBaselineTransform:
    """Verify the rigid_transform.txt was copied into the golden dir."""

    def test_rigid_transform_exists(self):
        assert (
            GOLDEN_DIR / "rigid_transform.txt"
        ).exists(), "rigid_transform.txt not found in legacy_baseline dir"

    def test_rigid_transform_is_4x4(self):
        txt_path = GOLDEN_DIR / "rigid_transform.txt"
        if not txt_path.exists():
            pytest.skip("rigid_transform.txt not yet copied")
        matrix = np.loadtxt(str(txt_path))
        assert matrix.shape == (4, 4), f"Expected 4x4 transform, got {matrix.shape}"
