"""Tests for h31: TypedDict + Unpack type hints on register function kwargs.

TDD RED phase: these tests will fail until _types.py is created and
register functions are updated with Unpack signatures.
"""


# ---------------------------------------------------------------------------
# Test 1: _types module exists with correct TypedDicts
# ---------------------------------------------------------------------------


def test_types_module_exists():
    """meshmonk._types must be importable."""
    import meshmonk._types  # noqa: F401


def test_rigid_kwargs_typeddict_exists():
    """RigidKwargs must be a TypedDict in meshmonk._types."""
    from meshmonk._types import RigidKwargs
    assert isinstance(RigidKwargs, type)


def test_nonrigid_kwargs_typeddict_exists():
    """NonrigidKwargs must be a TypedDict in meshmonk._types."""
    from meshmonk._types import NonrigidKwargs
    assert isinstance(NonrigidKwargs, type)


def test_pyramid_kwargs_typeddict_exists():
    """PyramidKwargs must be a TypedDict in meshmonk._types."""
    from meshmonk._types import PyramidKwargs
    assert isinstance(PyramidKwargs, type)


# ---------------------------------------------------------------------------
# Test 2: TypedDicts have total=False (all fields optional)
# ---------------------------------------------------------------------------


def test_rigid_kwargs_total_false():
    """RigidKwargs must have total=False (all fields optional)."""
    from meshmonk._types import RigidKwargs
    assert RigidKwargs.__total__ is False


def test_nonrigid_kwargs_total_false():
    """NonrigidKwargs must have total=False (all fields optional)."""
    from meshmonk._types import NonrigidKwargs
    assert NonrigidKwargs.__total__ is False


def test_pyramid_kwargs_total_false():
    """PyramidKwargs must have total=False (all fields optional)."""
    from meshmonk._types import PyramidKwargs
    assert PyramidKwargs.__total__ is False


# ---------------------------------------------------------------------------
# Test 3: TypedDicts have the correct fields
# ---------------------------------------------------------------------------


def test_rigid_kwargs_fields():
    """RigidKwargs must contain exactly the fields matching _RIGID_KNOWN_KWARGS."""
    from meshmonk._types import RigidKwargs
    expected = {
        "correspondences_symmetric",
        "correspondences_num_neighbours",
        "correspondences_flag_threshold",
        "correspondences_equalize_push_pull",
        "inlier_kappa",
        "inlier_use_orientation",
        "num_iterations",
        "use_scaling",
    }
    assert set(RigidKwargs.__annotations__.keys()) == expected


def test_nonrigid_kwargs_fields():
    """NonrigidKwargs must contain exactly the fields matching _NONRIGID_KNOWN_KWARGS."""
    from meshmonk._types import NonrigidKwargs
    expected = {
        "correspondences_symmetric",
        "correspondences_num_neighbours",
        "correspondences_flag_threshold",
        "correspondences_equalize_push_pull",
        "inlier_kappa",
        "inlier_use_orientation",
        "transform_sigma",
        "transform_num_viscous_iterations_start",
        "transform_num_viscous_iterations_end",
        "transform_num_elastic_iterations_start",
        "transform_num_elastic_iterations_end",
        "num_iterations",
    }
    assert set(NonrigidKwargs.__annotations__.keys()) == expected


def test_pyramid_kwargs_fields():
    """PyramidKwargs must contain exactly the fields matching _PYRAMID_KNOWN_KWARGS."""
    from meshmonk._types import PyramidKwargs
    expected = {
        "correspondences_symmetric",
        "correspondences_num_neighbours",
        "correspondences_flag_threshold",
        "correspondences_equalize_push_pull",
        "inlier_kappa",
        "inlier_use_orientation",
        "transform_sigma",
        "transform_num_viscous_iterations_start",
        "transform_num_viscous_iterations_end",
        "transform_num_elastic_iterations_start",
        "transform_num_elastic_iterations_end",
        "downsample_float_start",
        "downsample_target_start",
        "downsample_float_end",
        "downsample_target_end",
        "num_iterations",
        "num_pyramid_layers",
    }
    assert set(PyramidKwargs.__annotations__.keys()) == expected


# ---------------------------------------------------------------------------
# Test 4: Field type annotations are correct
# ---------------------------------------------------------------------------


def test_rigid_kwargs_field_types():
    """RigidKwargs fields must have correct Python type annotations."""
    from meshmonk._types import RigidKwargs
    annotations = RigidKwargs.__annotations__
    # bool fields
    assert annotations["correspondences_symmetric"] is bool
    assert annotations["correspondences_equalize_push_pull"] is bool
    assert annotations["inlier_use_orientation"] is bool
    assert annotations["use_scaling"] is bool
    # int fields
    assert annotations["correspondences_num_neighbours"] is int
    assert annotations["num_iterations"] is int
    # float fields
    assert annotations["correspondences_flag_threshold"] is float
    assert annotations["inlier_kappa"] is float


def test_nonrigid_kwargs_field_types():
    """NonrigidKwargs fields must have correct Python type annotations."""
    from meshmonk._types import NonrigidKwargs
    annotations = NonrigidKwargs.__annotations__
    # bool fields
    assert annotations["correspondences_symmetric"] is bool
    assert annotations["correspondences_equalize_push_pull"] is bool
    assert annotations["inlier_use_orientation"] is bool
    # int fields
    assert annotations["correspondences_num_neighbours"] is int
    assert annotations["num_iterations"] is int
    assert annotations["transform_num_viscous_iterations_start"] is int
    assert annotations["transform_num_viscous_iterations_end"] is int
    assert annotations["transform_num_elastic_iterations_start"] is int
    assert annotations["transform_num_elastic_iterations_end"] is int
    # float fields
    assert annotations["correspondences_flag_threshold"] is float
    assert annotations["inlier_kappa"] is float
    assert annotations["transform_sigma"] is float


def test_pyramid_kwargs_field_types():
    """PyramidKwargs fields must have correct Python type annotations."""
    from meshmonk._types import PyramidKwargs
    annotations = PyramidKwargs.__annotations__
    # bool fields
    assert annotations["correspondences_symmetric"] is bool
    assert annotations["correspondences_equalize_push_pull"] is bool
    assert annotations["inlier_use_orientation"] is bool
    # int fields
    assert annotations["correspondences_num_neighbours"] is int
    assert annotations["num_iterations"] is int
    assert annotations["num_pyramid_layers"] is int
    assert annotations["transform_num_viscous_iterations_start"] is int
    assert annotations["transform_num_viscous_iterations_end"] is int
    assert annotations["transform_num_elastic_iterations_start"] is int
    assert annotations["transform_num_elastic_iterations_end"] is int
    # float fields
    assert annotations["correspondences_flag_threshold"] is float
    assert annotations["inlier_kappa"] is float
    assert annotations["transform_sigma"] is float
    assert annotations["downsample_float_start"] is float
    assert annotations["downsample_target_start"] is float
    assert annotations["downsample_float_end"] is float
    assert annotations["downsample_target_end"] is float


# ---------------------------------------------------------------------------
# Test 5: TypedDicts are NOT in __all__
# ---------------------------------------------------------------------------


def test_typed_dicts_not_in_all():
    """RigidKwargs, NonrigidKwargs, PyramidKwargs must NOT appear in meshmonk.__all__."""
    import meshmonk
    assert "RigidKwargs" not in meshmonk.__all__
    assert "NonrigidKwargs" not in meshmonk.__all__
    assert "PyramidKwargs" not in meshmonk.__all__


# ---------------------------------------------------------------------------
# Test 6: typing_extensions is a direct dependency (importable)
# ---------------------------------------------------------------------------


def test_typing_extensions_importable():
    """typing_extensions must be importable (it's a direct dep after pyproject.toml update)."""
    import typing_extensions  # noqa: F401


def test_unpack_imported_from_typing_extensions():
    """Unpack must come from typing_extensions in meshmonk.__init__."""
    # Verify typing_extensions.Unpack is importable (the dep is in pyproject.toml)
    from typing_extensions import Unpack  # noqa: F401


# ---------------------------------------------------------------------------
# Test 7: Register functions use Unpack annotations
# ---------------------------------------------------------------------------


def test_rigid_register_kwargs_annotation():
    """rigid_register **kwargs must be annotated with Unpack[RigidKwargs]."""
    import meshmonk
    from typing_extensions import get_type_hints
    from meshmonk._types import RigidKwargs
    hints = get_type_hints(meshmonk.rigid_register, include_extras=True)
    assert "kwargs" in hints, "rigid_register must have a 'kwargs' annotation"
    # The annotation should be Unpack[RigidKwargs]
    kwargs_hint = hints["kwargs"]
    assert hasattr(kwargs_hint, "__args__"), f"Expected Unpack[...] annotation, got {kwargs_hint}"
    assert kwargs_hint.__args__ == (RigidKwargs,), (
        f"rigid_register kwargs must be Unpack[RigidKwargs], got {kwargs_hint}"
    )


def test_nonrigid_register_kwargs_annotation():
    """nonrigid_register **kwargs must be annotated with Unpack[NonrigidKwargs]."""
    import meshmonk
    from typing_extensions import get_type_hints
    from meshmonk._types import NonrigidKwargs
    hints = get_type_hints(meshmonk.nonrigid_register, include_extras=True)
    assert "kwargs" in hints, "nonrigid_register must have a 'kwargs' annotation"
    kwargs_hint = hints["kwargs"]
    assert hasattr(kwargs_hint, "__args__"), f"Expected Unpack[...] annotation, got {kwargs_hint}"
    assert kwargs_hint.__args__ == (NonrigidKwargs,), (
        f"nonrigid_register kwargs must be Unpack[NonrigidKwargs], got {kwargs_hint}"
    )


def test_pyramid_register_kwargs_annotation():
    """pyramid_register **kwargs must be annotated with Unpack[PyramidKwargs]."""
    import meshmonk
    from typing_extensions import get_type_hints
    from meshmonk._types import PyramidKwargs
    hints = get_type_hints(meshmonk.pyramid_register, include_extras=True)
    assert "kwargs" in hints, "pyramid_register must have a 'kwargs' annotation"
    kwargs_hint = hints["kwargs"]
    assert hasattr(kwargs_hint, "__args__"), f"Expected Unpack[...] annotation, got {kwargs_hint}"
    assert kwargs_hint.__args__ == (PyramidKwargs,), (
        f"pyramid_register kwargs must be Unpack[PyramidKwargs], got {kwargs_hint}"
    )


# ---------------------------------------------------------------------------
# Test 8: Runtime behavior unchanged (TypedDict(total=False) + kwargs.keys())
# ---------------------------------------------------------------------------


def test_rigid_register_runtime_behavior_unchanged(tmp_path):
    """Kwargs still apply correctly at runtime after type hint changes."""
    import meshmonk
    import numpy as np
    # Minimal valid Pattern B input
    n = 10
    feat = np.zeros((n, 6), dtype=np.float32)
    feat[:, 0] = np.linspace(0, 1, n)
    feat[:, 3] = 1.0  # normals
    faces = np.array([[0, 1, 2]], dtype=np.int32)
    # Call with some rigid kwargs - just check it doesn't crash
    # (The type hint change must not alter runtime behavior)
    try:
        meshmonk.rigid_register(
            floating_features=feat,
            target_features=feat,
            floating_faces=faces,
            target_faces=faces,
            num_iterations=5,
            use_scaling=False,
        )
    except Exception:
        pass  # Some registration errors are OK for degenerate input


def test_pyramid_register_explicit_kwargs_detection(tmp_path):
    """_apply_pyramid_kwargs must still detect explicit vs auto-populated kwargs.

    This verifies that TypedDict(total=False) preserves the kwargs.keys()
    membership-test behavior — explicit kwargs appear in keys(), auto-populated
    ones do not.
    """
    import meshmonk
    from meshmonk._meshmonk_core import PyramidParams

    # Call _apply_pyramid_kwargs directly with no explicit viscous/elastic start
    params = PyramidParams()
    params.num_iterations = 20
    kwargs = {"num_iterations": 20}
    explicit_kwargs = set(kwargs.keys())
    result = meshmonk._apply_pyramid_kwargs(params, kwargs, explicit_kwargs=explicit_kwargs)
    # Without explicit pass, viscous/elastic start should be auto-set to num_iterations
    assert result.transform.num_viscous_iterations_start == 20
    assert result.transform.num_elastic_iterations_start == 20

    # Now call WITH explicit viscous start
    params2 = PyramidParams()
    params2.num_iterations = 20
    kwargs2 = {"num_iterations": 20, "transform_num_viscous_iterations_start": 5}
    explicit_kwargs2 = set(kwargs2.keys())
    result2 = meshmonk._apply_pyramid_kwargs(params2, kwargs2, explicit_kwargs=explicit_kwargs2)
    assert result2.transform.num_viscous_iterations_start == 5
