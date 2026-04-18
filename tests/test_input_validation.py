"""Tests for Python input validation (bead workspace-hqz.1)."""
import numpy as np
import pytest
import meshmonk

# Small synthetic data for validation tests
def _make_features(n=10):
    """Create minimal (n,6) feature matrix and (n,3) faces."""
    rng = np.random.RandomState(42)
    verts = rng.randn(n, 3).astype(np.float32)
    normals = rng.randn(n, 3).astype(np.float32)
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


class TestEmptyFacesGuard:
    """Fix: rigid_registration must reject empty face arrays (DegenerateInput).

    nonrigid and pyramid already guard against empty faces; rigid did not.
    """

    def test_rigid_empty_floating_faces_raises(self):
        """rigid_register with 0-row floating_faces must raise MeshMonkError."""
        feat, _, flags = _make_features()
        empty_faces = np.zeros((0, 3), dtype=np.int32)
        real_faces = np.array([[0, 1, 2]], dtype=np.int32)
        with pytest.raises(meshmonk.MeshMonkError, match="DegenerateInput"):
            meshmonk.rigid_register(
                floating_features=feat,
                target_features=feat,
                floating_faces=empty_faces,
                target_faces=real_faces,
                floating_flags=flags,
                target_flags=flags,
            )

    def test_rigid_empty_target_faces_raises(self):
        """rigid_register with 0-row target_faces must raise MeshMonkError."""
        feat, _, flags = _make_features()
        real_faces = np.array([[0, 1, 2]], dtype=np.int32)
        empty_faces = np.zeros((0, 3), dtype=np.int32)
        with pytest.raises(meshmonk.MeshMonkError, match="DegenerateInput"):
            meshmonk.rigid_register(
                floating_features=feat,
                target_features=feat,
                floating_faces=real_faces,
                target_faces=empty_faces,
                floating_flags=flags,
                target_flags=flags,
            )


class TestPatternAFlagsNotSilentlyDropped:
    """Fix: passing floating_flags or target_flags with Pattern A must raise TypeError.

    Previously these were silently ignored (Pattern B detection didn't check flags).
    """

    def test_rigid_pattern_a_with_floating_flags_raises(self):
        """rigid_register: floating_flags with Pattern A should raise TypeError."""
        feat, faces, flags = _make_features()
        with pytest.raises(TypeError, match="Cannot mix Pattern A"):
            meshmonk.rigid_register(
                floating="some_path.obj",
                target="other_path.obj",
                floating_flags=flags,
            )

    def test_rigid_pattern_a_with_target_flags_raises(self):
        """rigid_register: target_flags with Pattern A should raise TypeError."""
        feat, faces, flags = _make_features()
        with pytest.raises(TypeError, match="Cannot mix Pattern A"):
            meshmonk.rigid_register(
                floating="some_path.obj",
                target="other_path.obj",
                target_flags=flags,
            )

    def test_nonrigid_pattern_a_with_floating_flags_raises(self):
        """nonrigid_register: floating_flags with Pattern A should raise TypeError."""
        feat, faces, flags = _make_features()
        with pytest.raises(TypeError, match="Cannot mix Pattern A"):
            meshmonk.nonrigid_register(
                floating="some_path.obj",
                target="other_path.obj",
                floating_flags=flags,
            )

    def test_pyramid_pattern_a_with_target_flags_raises(self):
        """pyramid_register: target_flags with Pattern A should raise TypeError."""
        feat, faces, flags = _make_features()
        with pytest.raises(TypeError, match="Cannot mix Pattern A"):
            meshmonk.pyramid_register(
                floating="some_path.obj",
                target="other_path.obj",
                target_flags=flags,
            )


class TestNonrigidAnnealingGuard:
    """Fix: NonrigidRegistration must guard annealing rate when viscous/elastic start=0.

    Without the guard: exp(log(end/0) / (iters-1)) = Inf/NaN — undefined behaviour.
    With the guard: rate = 1.0f (no annealing), matching PyramidNonrigidRegistration.

    We verify the fix via source-code inspection to directly capture the fix.
    """

    _SOURCE_FILE = "/workspace/library/src/NonrigidRegistration.cpp"

    def test_viscous_annealing_guard_present(self):
        """NonrigidRegistration.cpp must have a guard for zero viscous start value."""
        from pathlib import Path
        source = Path(self._SOURCE_FILE).read_text()
        # The fixed code must check _numViscousIterationsStart > 0 before computing rate
        assert "_numViscousIterationsStart > 0" in source, (
            "NonrigidRegistration.cpp is missing the guard for zero "
            "_numViscousIterationsStart. Without it, exp(log(end/0)/(iters-1)) "
            "produces Inf/NaN (undefined behaviour). Add the same guard that "
            "PyramidNonrigidRegistration.cpp has at lines 71-74."
        )

    def test_elastic_annealing_guard_present(self):
        """NonrigidRegistration.cpp must have a guard for zero elastic start value."""
        from pathlib import Path
        source = Path(self._SOURCE_FILE).read_text()
        # The fixed code must check _numElasticIterationsStart > 0 before computing rate
        assert "_numElasticIterationsStart > 0" in source, (
            "NonrigidRegistration.cpp is missing the guard for zero "
            "_numElasticIterationsStart. Without it, exp(log(end/0)/(iters-1)) "
            "produces Inf/NaN (undefined behaviour). Add the same guard that "
            "PyramidNonrigidRegistration.cpp has at lines 75-78."
        )
