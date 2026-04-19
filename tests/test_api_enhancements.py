"""Tests for API enhancements (bead workspace-hqz.6)."""

import os
import subprocess
import sys

import numpy as np
import pytest


class TestMainModule:
    def test_python_m_meshmonk_help(self):
        """python -m meshmonk --help should print CLI help."""
        env = dict(os.environ)  # inherit current environment (includes LD_LIBRARY_PATH)
        result = subprocess.run(
            [sys.executable, "-m", "meshmonk", "--help"],
            capture_output=True,
            text=True,
            timeout=30,
            env=env,
        )
        assert result.returncode == 0
        assert "meshmonk" in result.stdout.lower()


class TestReadOnlyTransformMatrix:
    def test_transform_matrix_not_writable(self):
        """RigidTransform.matrix should be read-only from Python."""
        from meshmonk._meshmonk_core import RigidTransform

        t = RigidTransform()
        # Reading should work
        m = t.matrix
        assert m.shape == (4, 4)
        # Writing should raise
        with pytest.raises((AttributeError, TypeError)):
            t.matrix = np.eye(4, dtype=np.float32)


class TestPatternAFlags:
    def test_mesh_flags_respected(self):
        """Pattern A should use mesh.flags if available."""
        pytest.importorskip("trimesh")
        import trimesh

        # Create a simple mesh with custom flags
        vertices = np.array(
            [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64
        )
        faces = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]], dtype=np.int64)
        mesh = trimesh.Trimesh(vertices=vertices, faces=faces)

        # Add flags attribute
        mesh.flags = np.array([1.0, 1.0, 0.0, 1.0], dtype=np.float32)

        # Call _mesh_to_arrays and check flags are used
        from meshmonk import _mesh_to_arrays

        features, f, flags = _mesh_to_arrays(mesh)
        np.testing.assert_array_equal(flags, mesh.flags)

    def test_mesh_without_flags_gets_ones(self):
        """Pattern A without .flags should get all-ones flags."""
        pytest.importorskip("trimesh")
        import trimesh

        vertices = np.array(
            [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64
        )
        faces = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]], dtype=np.int64)
        mesh = trimesh.Trimesh(vertices=vertices, faces=faces)

        from meshmonk import _mesh_to_arrays

        features, f, flags = _mesh_to_arrays(mesh)
        np.testing.assert_array_equal(flags, np.ones(4, dtype=np.float32))

    def test_mesh_flags_wrong_shape_raises(self):
        """mesh.flags with wrong shape should raise ValueError."""
        pytest.importorskip("trimesh")
        import trimesh

        vertices = np.array(
            [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64
        )
        faces = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]], dtype=np.int64)
        mesh = trimesh.Trimesh(vertices=vertices, faces=faces)
        mesh.flags = np.array([1.0, 1.0], dtype=np.float32)  # wrong shape

        from meshmonk import _mesh_to_arrays

        with pytest.raises(ValueError, match="flags shape"):
            _mesh_to_arrays(mesh)


class TestPrepareArrays:
    def _make_tetrahedron(self):
        import trimesh

        verts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)
        faces = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]], dtype=np.int32)
        return trimesh.Trimesh(vertices=verts, faces=faces, process=False)

    def test_pattern_a_returns_correct_shape(self):
        from meshmonk import _prepare_arrays

        pytest.importorskip("trimesh")
        mesh = self._make_tetrahedron()
        feat_f, feat_t, faces_f, faces_t, flags_f, flags_t = _prepare_arrays(mesh, mesh)
        assert feat_f.shape == (4, 6)
        assert feat_f.dtype == np.float32
        assert faces_f.shape == (4, 3)
        assert flags_f.shape == (4,)
        np.testing.assert_array_equal(flags_f, np.ones(4, dtype=np.float32))

    def test_pattern_a_flags_override(self):
        from meshmonk import _prepare_arrays

        pytest.importorskip("trimesh")
        mesh = self._make_tetrahedron()
        custom = np.array([1.0, 0.5, 0.0, 1.0], dtype=np.float32)
        _, _, _, _, flags_f, _ = _prepare_arrays(mesh, mesh, floating_flags=custom)
        np.testing.assert_array_equal(flags_f, custom)

    def test_pattern_b_returns_correct_shape(self):
        from meshmonk import _prepare_arrays

        feat = np.random.randn(4, 6).astype(np.float32)
        faces = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]], dtype=np.int32)
        feat_f, feat_t, faces_f, faces_t, flags_f, flags_t = _prepare_arrays(
            None,
            None,
            floating_features=feat,
            target_features=feat,
            floating_faces=faces,
            target_faces=faces,
        )
        assert feat_f.shape == (4, 6)
        assert flags_f.shape == (4,)
        np.testing.assert_array_equal(flags_f, np.ones(4, dtype=np.float32))

    def test_pattern_b_missing_required_raises(self):
        from meshmonk import _prepare_arrays

        feat = np.random.randn(4, 6).astype(np.float32)
        with pytest.raises(ValueError, match="Missing"):
            _prepare_arrays(None, None, floating_features=feat)


class TestRigidParamsKwarg:
    def _make_mesh(self):
        """Return a mesh large enough for rigid pre-alignment to work reliably."""
        import trimesh

        return trimesh.creation.icosphere(subdivisions=1)

    def test_nonrigid_rigid_params_empty_dict_runs(self):
        import meshmonk

        pytest.importorskip("trimesh")
        mesh = self._make_mesh()
        result = meshmonk.nonrigid_register(
            floating=mesh, target=mesh, rigid_params={}, num_iterations=3
        )
        assert result.aligned_features.shape == (mesh.vertices.shape[0], 6)

    def test_pyramid_rigid_params_empty_dict_runs(self):
        import meshmonk

        pytest.importorskip("trimesh")
        mesh = self._make_mesh()
        result = meshmonk.pyramid_register(
            floating=mesh,
            target=mesh,
            rigid_params={},
            num_iterations=4,
            num_pyramid_layers=1,
        )
        assert result.aligned_features.shape == (mesh.vertices.shape[0], 6)

    def test_rigid_params_none_preserves_behavior(self):
        import meshmonk

        pytest.importorskip("trimesh")
        mesh = self._make_mesh()
        result_default = meshmonk.nonrigid_register(
            floating=mesh, target=mesh, num_iterations=3
        )
        result_explicit_none = meshmonk.nonrigid_register(
            floating=mesh, target=mesh, rigid_params=None, num_iterations=3
        )
        assert (
            result_default.aligned_features.shape
            == result_explicit_none.aligned_features.shape
        )
