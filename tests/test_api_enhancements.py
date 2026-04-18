"""Tests for v0.2 API enhancements (bead workspace-hqz.6)."""
import subprocess
import sys

import numpy as np
import pytest


class TestMainModule:
    def test_python_m_meshmonk_help(self):
        """python -m meshmonk --help should print CLI help."""
        result = subprocess.run(
            [sys.executable, "-m", "meshmonk", "--help"],
            capture_output=True, text=True, timeout=30,
            env={**__import__('os').environ, "LD_LIBRARY_PATH": "/workspace/.venv/lib/python3.11/site-packages/lib"},
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
        vertices = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]], dtype=np.float64)
        faces = np.array([[0,1,2],[0,1,3],[0,2,3],[1,2,3]], dtype=np.int64)
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

        vertices = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]], dtype=np.float64)
        faces = np.array([[0,1,2],[0,1,3],[0,2,3],[1,2,3]], dtype=np.int64)
        mesh = trimesh.Trimesh(vertices=vertices, faces=faces)

        from meshmonk import _mesh_to_arrays
        features, f, flags = _mesh_to_arrays(mesh)
        np.testing.assert_array_equal(flags, np.ones(4, dtype=np.float32))

    def test_mesh_flags_wrong_shape_raises(self):
        """mesh.flags with wrong shape should raise ValueError."""
        pytest.importorskip("trimesh")
        import trimesh

        vertices = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]], dtype=np.float64)
        faces = np.array([[0,1,2],[0,1,3],[0,2,3],[1,2,3]], dtype=np.int64)
        mesh = trimesh.Trimesh(vertices=vertices, faces=faces)
        mesh.flags = np.array([1.0, 1.0], dtype=np.float32)  # wrong shape

        from meshmonk import _mesh_to_arrays
        with pytest.raises(ValueError, match="flags shape"):
            _mesh_to_arrays(mesh)
