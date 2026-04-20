"""Tests for kwargs deduplication refactoring.

These tests verify that the refactored _apply_*_kwargs functions produce
identical behavior to the pre-refactoring versions. They test:
1. Shared helpers exist and are callable
2. All kwargs still apply correctly to each registration type
3. Unknown kwargs still raise TypeError
4. End-to-end registration still works with kwargs
"""

import pytest
import meshmonk
from meshmonk._meshmonk_core import (
    RigidParams,
    NonrigidParams,
    PyramidParams,
)


# ---------------------------------------------------------------------------
# Verify shared helpers exist
# ---------------------------------------------------------------------------


def test_apply_shared_kwargs_exists():
    """_apply_shared_kwargs must exist as an internal helper."""
    assert hasattr(
        meshmonk, "_apply_shared_kwargs"
    ), "Expected meshmonk._apply_shared_kwargs to exist after refactoring"
    assert callable(meshmonk._apply_shared_kwargs)


def test_apply_transform_kwargs_exists():
    """_apply_transform_kwargs must exist as an internal helper."""
    assert hasattr(
        meshmonk, "_apply_transform_kwargs"
    ), "Expected meshmonk._apply_transform_kwargs to exist after refactoring"
    assert callable(meshmonk._apply_transform_kwargs)


# ---------------------------------------------------------------------------
# Shared kwargs apply correctly to all param types
# ---------------------------------------------------------------------------


class TestSharedKwargsApplyToAllTypes:
    """Correspondence + inlier kwargs must apply identically to rigid, nonrigid, pyramid."""

    SHARED_KWARGS = {
        "correspondences_symmetric": False,
        "correspondences_num_neighbours": 5,
        "correspondences_flag_threshold": 0.8,
        "correspondences_equalize_push_pull": True,
        "inlier_kappa": 8.0,
        "inlier_use_orientation": False,
    }

    def test_rigid_accepts_all_shared_kwargs(self):
        params = RigidParams()
        params = meshmonk._apply_rigid_kwargs(params, self.SHARED_KWARGS.copy())
        assert not params.correspondences.symmetric
        assert params.correspondences.num_neighbours == 5
        assert params.correspondences.flag_threshold == pytest.approx(0.8)
        assert params.correspondences.equalize_push_pull
        assert params.inliers.kappa == pytest.approx(8.0)
        assert not params.inliers.use_orientation

    def test_nonrigid_accepts_all_shared_kwargs(self):
        params = NonrigidParams()
        params = meshmonk._apply_nonrigid_kwargs(params, self.SHARED_KWARGS.copy())
        assert not params.correspondences.symmetric
        assert params.correspondences.num_neighbours == 5
        assert params.correspondences.flag_threshold == pytest.approx(0.8)
        assert params.correspondences.equalize_push_pull
        assert params.inliers.kappa == pytest.approx(8.0)
        assert not params.inliers.use_orientation

    def test_pyramid_accepts_all_shared_kwargs(self):
        params = PyramidParams()
        params = meshmonk._apply_pyramid_kwargs(
            params,
            self.SHARED_KWARGS.copy(),
            explicit_kwargs=set(self.SHARED_KWARGS.keys()),
        )
        assert not params.correspondences.symmetric
        assert params.correspondences.num_neighbours == 5
        assert params.correspondences.flag_threshold == pytest.approx(0.8)
        assert params.correspondences.equalize_push_pull
        assert params.inliers.kappa == pytest.approx(8.0)
        assert not params.inliers.use_orientation


# ---------------------------------------------------------------------------
# Transform kwargs shared by nonrigid and pyramid
# ---------------------------------------------------------------------------


class TestTransformKwargs:
    """Transform kwargs must apply to both nonrigid and pyramid."""

    TRANSFORM_KWARGS = {
        "transform_sigma": 5.0,
        "transform_num_viscous_iterations_start": 100,
        "transform_num_viscous_iterations_end": 10,
        "transform_num_elastic_iterations_start": 100,
        "transform_num_elastic_iterations_end": 10,
    }

    def test_nonrigid_accepts_transform_kwargs(self):
        params = NonrigidParams()
        params = meshmonk._apply_nonrigid_kwargs(params, self.TRANSFORM_KWARGS.copy())
        assert params.transform.sigma == pytest.approx(5.0)
        assert params.transform.num_viscous_iterations_start == 100
        assert params.transform.num_viscous_iterations_end == 10
        assert params.transform.num_elastic_iterations_start == 100
        assert params.transform.num_elastic_iterations_end == 10

    def test_pyramid_accepts_transform_kwargs(self):
        params = PyramidParams()
        params = meshmonk._apply_pyramid_kwargs(
            params,
            self.TRANSFORM_KWARGS.copy(),
            explicit_kwargs=set(self.TRANSFORM_KWARGS.keys()),
        )
        assert params.transform.sigma == pytest.approx(5.0)
        assert params.transform.num_viscous_iterations_start == 100
        assert params.transform.num_viscous_iterations_end == 10
        assert params.transform.num_elastic_iterations_start == 100
        assert params.transform.num_elastic_iterations_end == 10


# ---------------------------------------------------------------------------
# Type-specific kwargs still work
# ---------------------------------------------------------------------------


def test_rigid_specific_kwargs():
    params = RigidParams()
    params = meshmonk._apply_rigid_kwargs(
        params,
        {
            "num_iterations": 20,
            "use_scaling": True,
        },
    )
    assert params.num_iterations == 20
    assert params.use_scaling


def test_nonrigid_specific_kwargs():
    params = NonrigidParams()
    params = meshmonk._apply_nonrigid_kwargs(params, {"num_iterations": 15})
    assert params.num_iterations == 15


def test_pyramid_specific_kwargs():
    params = PyramidParams()
    params = meshmonk._apply_pyramid_kwargs(
        params,
        {
            "downsample_float_start": 500,
            "downsample_target_start": 600,
            "downsample_float_end": 100,
            "downsample_target_end": 200,
            "num_iterations": 10,
            "num_pyramid_layers": 3,
        },
        explicit_kwargs={
            "downsample_float_start",
            "downsample_target_start",
            "downsample_float_end",
            "downsample_target_end",
            "num_iterations",
            "num_pyramid_layers",
        },
    )
    assert params.downsample.float_start == 500
    assert params.downsample.target_start == 600
    assert params.downsample.float_end == 100
    assert params.downsample.target_end == 200
    assert params.num_iterations == 10
    assert params.num_pyramid_layers == 3


# ---------------------------------------------------------------------------
# MATLAB convention: viscous/elastic start auto-populated
# ---------------------------------------------------------------------------


def test_pyramid_matlab_convention_auto_viscous_elastic():
    """When user sets num_iterations but NOT viscous/elastic start,
    those should auto-populate to num_iterations value."""
    params = PyramidParams()
    params = meshmonk._apply_pyramid_kwargs(
        params,
        {"num_iterations": 42},
        explicit_kwargs={"num_iterations"},
    )
    assert params.num_iterations == 42
    assert params.transform.num_viscous_iterations_start == 42
    assert params.transform.num_elastic_iterations_start == 42


def test_pyramid_matlab_convention_respects_explicit():
    """When user explicitly sets viscous_start, it should NOT be overridden."""
    params = PyramidParams()
    params = meshmonk._apply_pyramid_kwargs(
        params,
        {"num_iterations": 42, "transform_num_viscous_iterations_start": 99},
        explicit_kwargs={"num_iterations", "transform_num_viscous_iterations_start"},
    )
    assert params.num_iterations == 42
    assert (
        params.transform.num_viscous_iterations_start == 99
    )  # explicit, not overridden
    assert params.transform.num_elastic_iterations_start == 42  # auto


# ---------------------------------------------------------------------------
# Unknown kwargs still raise TypeError
# ---------------------------------------------------------------------------


def test_rigid_unknown_kwarg_raises():
    params = RigidParams()
    with pytest.raises(TypeError, match="unexpected keyword"):
        meshmonk._apply_rigid_kwargs(params, {"bogus_param": 1})


def test_nonrigid_unknown_kwarg_raises():
    params = NonrigidParams()
    with pytest.raises(TypeError, match="unexpected keyword"):
        meshmonk._apply_nonrigid_kwargs(params, {"bogus_param": 1})


def test_pyramid_unknown_kwarg_raises():
    params = PyramidParams()
    with pytest.raises(TypeError, match="unexpected keyword"):
        meshmonk._apply_pyramid_kwargs(
            params, {"bogus_param": 1}, explicit_kwargs=set()
        )


# ---------------------------------------------------------------------------
# End-to-end: registration with kwargs still works
# ---------------------------------------------------------------------------


def test_rigid_register_with_kwargs_e2e():
    pytest.importorskip("trimesh")
    import trimesh

    m = trimesh.creation.icosphere(subdivisions=0)
    result = meshmonk.rigid_register(
        floating=m,
        target=m,
        num_iterations=3,
        correspondences_symmetric=True,
        inlier_kappa=10.0,
    )
    assert result.aligned_features.shape == (len(m.vertices), 6)


def test_nonrigid_register_with_kwargs_e2e():
    pytest.importorskip("trimesh")
    import trimesh

    m = trimesh.creation.icosphere(subdivisions=0)
    result = meshmonk.nonrigid_register(
        floating=m,
        target=m,
        num_iterations=3,
        correspondences_num_neighbours=4,
        transform_sigma=4.0,
    )
    assert result.aligned_features.shape == (len(m.vertices), 6)


def test_pyramid_register_with_kwargs_e2e():
    pytest.importorskip("trimesh")
    import trimesh

    m = trimesh.creation.icosphere(subdivisions=0)
    result = meshmonk.pyramid_register(
        floating=m,
        target=m,
        num_iterations=4,
        num_pyramid_layers=1,
        correspondences_flag_threshold=0.95,
        inlier_use_orientation=False,
    )
    assert result.aligned_features.shape == (len(m.vertices), 6)
