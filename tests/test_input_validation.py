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

    def test_pattern_b_none_flags_defaults_to_ones(self):
        feat, faces, _ = _make_features()
        result = meshmonk.rigid_register(
            floating_features=feat,
            target_features=feat,
            floating_faces=faces,
            target_faces=faces,
            floating_flags=None,
            target_flags=None,
        )
        assert result is not None


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

    _SOURCE_FILE = str(
        __import__("pathlib").Path(__file__).resolve().parent.parent
        / "library"
        / "src"
        / "NonrigidRegistration.cpp"
    )

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


# ---------------------------------------------------------------------------
# New edge-case tests (bead meshmonk-modernization-ald)
# ---------------------------------------------------------------------------


class TestZeroIterations:
    """Document behavior of num_iterations=0 across registration functions.

    Current behavior (documented, not fixed here):
    - rigid_register: returns unchanged input (iterations_run=0)
    - nonrigid_register: raises MeshMonkError(DegenerateInput)
    - pyramid_register: raises MeshMonkError(DegenerateInput)

    The divergence between rigid and nonrigid/pyramid is surprising;
    see bead meshmonk-modernization-ald for context.
    """

    def test_rigid_zero_iterations_returns_unchanged_input(self):
        """rigid_register with num_iterations=0 returns input unchanged (iterations_run=0)."""
        feat, faces, flags = _make_features()
        result = meshmonk.rigid_register(
            floating_features=feat,
            target_features=feat,
            floating_faces=faces,
            target_faces=faces,
            floating_flags=flags,
            target_flags=flags,
            num_iterations=0,
        )
        assert result.iterations_run == 0
        np.testing.assert_array_equal(result.aligned_features, feat)

    @pytest.mark.parametrize(
        "register_fn, fn_name",
        [
            (meshmonk.nonrigid_register, "nonrigid_register"),
            (meshmonk.pyramid_register, "pyramid_register"),
        ],
    )
    def test_nonrigid_pyramid_zero_iterations_raises_degenerate_input(
        self, register_fn, fn_name
    ):
        """nonrigid/pyramid_register with num_iterations=0 raises MeshMonkError(DegenerateInput).

        Current behavior: the C++ layer raises DegenerateInput when num_iterations=0
        for nonrigid and pyramid — unlike rigid, which silently returns unchanged input.
        """
        feat, faces, flags = _make_features()
        with pytest.raises(meshmonk.MeshMonkError, match="DegenerateInput"):
            register_fn(
                floating_features=feat,
                target_features=feat,
                floating_faces=faces,
                target_faces=faces,
                floating_flags=flags,
                target_flags=flags,
                num_iterations=0,
            )


class TestDegenerateMeshes:
    """Document behavior of degenerate mesh inputs.

    All tested via rigid_register (the lightest pipeline).
    """

    def test_single_vertex_raises_degenerate_input(self):
        """A single-vertex floating mesh raises MeshMonkError(DegenerateInput)."""
        feat, faces, flags = _make_features()
        single_vert = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, 1.0]], dtype=np.float32)
        single_flags = np.ones(1, dtype=np.float32)
        empty_faces = np.zeros((0, 3), dtype=np.int32)
        with pytest.raises(meshmonk.MeshMonkError, match="DegenerateInput"):
            meshmonk.rigid_register(
                floating_features=single_vert,
                target_features=feat,
                floating_faces=empty_faces,
                target_faces=faces,
                floating_flags=single_flags,
                target_flags=flags,
            )

    def test_empty_vertex_raises_degenerate_input(self):
        """A zero-row floating mesh raises MeshMonkError(DegenerateInput)."""
        feat, faces, flags = _make_features()
        empty_feat = np.zeros((0, 6), dtype=np.float32)
        empty_flags = np.zeros(0, dtype=np.float32)
        empty_faces = np.zeros((0, 3), dtype=np.int32)
        with pytest.raises(meshmonk.MeshMonkError, match="DegenerateInput"):
            meshmonk.rigid_register(
                floating_features=empty_feat,
                target_features=feat,
                floating_faces=empty_faces,
                target_faces=faces,
                floating_flags=empty_flags,
                target_flags=flags,
            )

    def test_single_face_mesh_raises_insufficient_inliers(self):
        """A 3-vertex / 1-face floating mesh raises MeshMonkError(InsufficientInliers).

        The C++ ICP correspondences require >=4 inlier vertices; 3 is too few.
        """
        feat, faces, flags = _make_features()
        rng = np.random.RandomState(42)
        sv = rng.randn(3, 3).astype(np.float32)
        sn = rng.randn(3, 3).astype(np.float32)
        # Normalize normals so they are non-degenerate
        sn /= np.linalg.norm(sn, axis=1, keepdims=True)
        sf = np.hstack([sv, sn])
        sfac = np.array([[0, 1, 2]], dtype=np.int32)
        sfl = np.ones(3, dtype=np.float32)
        with pytest.raises(meshmonk.MeshMonkError, match="InsufficientInliers"):
            meshmonk.rigid_register(
                floating_features=sf,
                target_features=feat,
                floating_faces=sfac,
                target_faces=faces,
                floating_flags=sfl,
                target_flags=flags,
            )


class TestNonFiniteNormals:
    """Document behavior when the features matrix contains NaN or inf normals.

    Current behavior: the C++ layer eventually raises MeshMonkError(DecompositionFailed)
    after running many failed eigenvector decompositions — no early validation at the
    Python boundary.

    FIXME(meshmonk-modernization-2ke): non-finite normals should be caught at the Python
    boundary in _prepare_arrays and raise a clear ValueError instead of letting the C++
    layer spam warning lines and raise DecompositionFailed.
    """

    @pytest.mark.parametrize(
        "nonfinite_value, label",
        [
            (np.nan, "NaN"),
            (np.inf, "inf"),
        ],
    )
    def test_nonfinite_normals_raises_decomposition_failed(
        self, nonfinite_value, label
    ):
        """NaN or inf normals in floating_features raise MeshMonkError(DecompositionFailed).

        FIXME(meshmonk-modernization-2ke): current behavior — the C++ layer raises
        DecompositionFailed after many failed eigenvector decompositions instead of
        the Python boundary validating and raising ValueError early.
        """
        feat, faces, flags = _make_features()
        bad_feat = feat.copy()
        bad_feat[:, 3:] = nonfinite_value  # corrupt normals columns only
        with pytest.raises(meshmonk.MeshMonkError, match="DecompositionFailed"):
            meshmonk.rigid_register(
                floating_features=bad_feat,
                target_features=feat,
                floating_faces=faces,
                target_faces=faces,
                floating_flags=flags,
                target_flags=flags,
            )


class TestDtypeCoercion:
    """Document silent dtype coercion in _prepare_arrays (Pattern B).

    Current behavior: float64 features/flags and int64 faces are silently coerced
    to float32/int32 via np.asarray(..., dtype=...) with no warning emitted.
    The output is always float32 regardless of input dtype.
    """

    @pytest.mark.parametrize(
        "feat_dtype, faces_dtype, flags_dtype",
        [
            (np.float64, np.int32, np.float32),   # float64 features only
            (np.float32, np.int64, np.float32),   # int64 faces only
            (np.float32, np.int32, np.float64),   # float64 flags only
            (np.float64, np.int64, np.float64),   # all wider dtypes
        ],
    )
    def test_dtype_coercion_succeeds_silently(
        self, feat_dtype, faces_dtype, flags_dtype
    ):
        """Wider dtypes are silently coerced to float32/int32 without warning or error.

        Current behavior: np.asarray(x, dtype=float32/int32) coerces silently —
        no DeprecationWarning, no TypeError raised.
        """
        import warnings

        feat, faces, flags = _make_features()
        feat_in = feat.astype(feat_dtype)
        faces_in = faces.astype(faces_dtype)
        flags_in = flags.astype(flags_dtype)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = meshmonk.rigid_register(
                floating_features=feat_in,
                target_features=feat_in,
                floating_faces=faces_in,
                target_faces=faces_in,
                floating_flags=flags_in,
                target_flags=flags_in,
            )

        # No warnings should be emitted for dtype coercion
        dtype_warnings = [x for x in w if "dtype" in str(x.message).lower()]
        assert dtype_warnings == [], f"Unexpected dtype warnings: {dtype_warnings}"

        # Output is always float32 regardless of input dtype
        assert result.aligned_features.dtype == np.float32
