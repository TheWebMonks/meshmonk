"""Tests for Python input validation (bead workspace-hqz.1)."""
import numpy as np
import pytest
import meshmonk

# Small synthetic data for validation tests
def _make_features(n=10):
    """Create minimal (n,6) feature matrix and (n,3) faces."""
    verts = np.random.randn(n, 3).astype(np.float32)
    normals = np.random.randn(n, 3).astype(np.float32)
    features = np.hstack([verts, normals])
    # Minimal triangulation
    faces = np.array([[0, 1, 2]], dtype=np.int32)
    flags = np.ones(n, dtype=np.float32)
    return features, faces, flags


class TestPatternBNoneGuard:
    def test_rigid_none_floating_features_raises(self):
        feat, faces, flags = _make_features()
        with pytest.raises(ValueError, match="Pattern B requires"):
            meshmonk.rigid_register(
                floating_features=None,
                target_features=feat,
                floating_faces=faces,
                target_faces=faces,
                floating_flags=flags,
                target_flags=flags,
            )

    def test_nonrigid_none_target_faces_raises(self):
        feat, faces, flags = _make_features()
        with pytest.raises(ValueError, match="Pattern B requires"):
            meshmonk.nonrigid_register(
                floating_features=feat,
                target_features=feat,
                floating_faces=faces,
                target_faces=None,
                floating_flags=flags,
                target_flags=flags,
            )

    def test_pyramid_none_flags_raises(self):
        feat, faces, flags = _make_features()
        with pytest.raises(ValueError, match="Pattern B requires"):
            meshmonk.pyramid_register(
                floating_features=feat,
                target_features=feat,
                floating_faces=faces,
                target_faces=faces,
                floating_flags=None,
                target_flags=flags,
            )


class TestMixedPatternRejection:
    def test_rigid_mixed_raises(self):
        feat, faces, flags = _make_features()
        with pytest.raises(TypeError, match="Cannot mix Pattern A"):
            meshmonk.rigid_register(
                floating="some_path.obj",
                target="other_path.obj",
                floating_features=feat,
            )

    def test_nonrigid_mixed_raises(self):
        feat, faces, flags = _make_features()
        with pytest.raises(TypeError, match="Cannot mix Pattern A"):
            meshmonk.nonrigid_register(
                floating="some_path.obj",
                target="other_path.obj",
                target_features=feat,
            )


class TestNormalsOverrideValidation:
    def test_wrong_shape_raises(self):
        """normals_override with wrong shape should raise ValueError."""
        pytest.importorskip("trimesh")
        from pathlib import Path
        template = Path("/workspace/data/Template.obj")
        if not template.exists():
            pytest.skip("Template.obj not available")
        bad_normals = np.random.randn(5, 3).astype(np.float32)  # wrong N
        with pytest.raises(ValueError, match="normals_override shape"):
            meshmonk.rigid_register(
                floating=str(template),
                target=str(template),
                normals=bad_normals,
            )
