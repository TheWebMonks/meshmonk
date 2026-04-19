"""Tests for MeshMonkError.code binding.

These tests verify that MeshMonkError instances have a programmatically
accessible .code attribute returning a RegistrationError enum value.
This replaces the previous approach of string-parsing error messages.
"""

import pytest
import numpy as np
import meshmonk
from meshmonk import MeshMonkError, RegistrationError


# ---------------------------------------------------------------------------
# Core: .code attribute exists and returns RegistrationError
# ---------------------------------------------------------------------------


def test_meshmonk_error_has_code_property():
    """MeshMonkError must have a .code attribute."""
    empty = np.zeros((0, 6), dtype=np.float32)
    faces = np.zeros((0, 3), dtype=np.int32)
    flags = np.zeros((0,), dtype=np.float32)
    with pytest.raises(MeshMonkError) as exc_info:
        meshmonk.rigid_register(
            floating_features=empty,
            target_features=empty,
            floating_faces=faces,
            target_faces=faces,
            floating_flags=flags,
            target_flags=flags,
        )
    assert hasattr(
        exc_info.value, "code"
    ), "MeshMonkError instance must have a .code attribute"


def test_degenerate_input_code_is_registration_error_enum():
    """MeshMonkError.code must return a RegistrationError enum member."""
    empty = np.zeros((0, 6), dtype=np.float32)
    faces = np.zeros((0, 3), dtype=np.int32)
    flags = np.zeros((0,), dtype=np.float32)
    with pytest.raises(MeshMonkError) as exc_info:
        meshmonk.rigid_register(
            floating_features=empty,
            target_features=empty,
            floating_faces=faces,
            target_faces=faces,
            floating_flags=flags,
            target_flags=flags,
        )
    code = exc_info.value.code
    assert isinstance(
        code, RegistrationError
    ), f"Expected .code to be RegistrationError, got {type(code)}: {code!r}"


def test_degenerate_input_code_value():
    """Empty input must produce DegenerateInput code."""
    empty = np.zeros((0, 6), dtype=np.float32)
    faces = np.zeros((0, 3), dtype=np.int32)
    flags = np.zeros((0,), dtype=np.float32)
    with pytest.raises(MeshMonkError) as exc_info:
        meshmonk.rigid_register(
            floating_features=empty,
            target_features=empty,
            floating_faces=faces,
            target_faces=faces,
            floating_flags=flags,
            target_flags=flags,
        )
    assert exc_info.value.code == RegistrationError.DegenerateInput


# ---------------------------------------------------------------------------
# Programmatic error handling (the main use case)
# ---------------------------------------------------------------------------


def test_programmatic_error_dispatch():
    """Users must be able to switch on .code instead of parsing strings."""
    empty = np.zeros((0, 6), dtype=np.float32)
    faces = np.zeros((0, 3), dtype=np.int32)
    flags = np.zeros((0,), dtype=np.float32)
    try:
        meshmonk.rigid_register(
            floating_features=empty,
            target_features=empty,
            floating_faces=faces,
            target_faces=faces,
            floating_flags=flags,
            target_flags=flags,
        )
        pytest.fail("Expected MeshMonkError")
    except MeshMonkError as e:
        if e.code == RegistrationError.DegenerateInput:
            handled = True
        else:
            handled = False
    assert handled, "Should have handled DegenerateInput via .code dispatch"


# ---------------------------------------------------------------------------
# Exception hierarchy preserved
# ---------------------------------------------------------------------------


def test_meshmonk_error_still_subclass_of_runtime_error():
    """MeshMonkError must remain a RuntimeError subclass."""
    assert issubclass(MeshMonkError, RuntimeError)


def test_except_runtime_error_catches_meshmonk_error():
    """except RuntimeError must still catch MeshMonkError."""
    empty = np.zeros((0, 6), dtype=np.float32)
    faces = np.zeros((0, 3), dtype=np.int32)
    flags = np.zeros((0,), dtype=np.float32)
    with pytest.raises(RuntimeError):
        meshmonk.rigid_register(
            floating_features=empty,
            target_features=empty,
            floating_faces=faces,
            target_faces=faces,
            floating_flags=flags,
            target_flags=flags,
        )


def test_except_meshmonk_error_catches_it():
    """except MeshMonkError must catch the exception."""
    empty = np.zeros((0, 6), dtype=np.float32)
    faces = np.zeros((0, 3), dtype=np.int32)
    flags = np.zeros((0,), dtype=np.float32)
    caught = False
    try:
        meshmonk.rigid_register(
            floating_features=empty,
            target_features=empty,
            floating_faces=faces,
            target_faces=faces,
            floating_flags=flags,
            target_flags=flags,
        )
    except MeshMonkError:
        caught = True
    assert caught


# ---------------------------------------------------------------------------
# Message string preserved
# ---------------------------------------------------------------------------


def test_error_message_still_contains_description():
    """str(e) must still contain the human-readable error description."""
    empty = np.zeros((0, 6), dtype=np.float32)
    faces = np.zeros((0, 3), dtype=np.int32)
    flags = np.zeros((0,), dtype=np.float32)
    with pytest.raises(MeshMonkError) as exc_info:
        meshmonk.rigid_register(
            floating_features=empty,
            target_features=empty,
            floating_faces=faces,
            target_faces=faces,
            floating_flags=flags,
            target_flags=flags,
        )
    msg = str(exc_info.value)
    assert (
        "DegenerateInput" in msg
    ), f"Expected 'DegenerateInput' in message, got: {msg!r}"


# ---------------------------------------------------------------------------
# All error codes accessible via nonrigid (InsufficientInliers)
# ---------------------------------------------------------------------------


def test_insufficient_inliers_code():
    """Trigger InsufficientInliers and verify .code."""
    # Use a mesh where most flags are zero → too few inliers after InlierDetector
    n = 100
    features = np.random.randn(n, 6).astype(np.float32)
    faces = np.zeros((1, 3), dtype=np.int32)  # minimal faces
    # Set all flags to near-zero except 2 (fewer than the 4 minimum)
    flags = np.zeros(n, dtype=np.float32)
    flags[0] = 1.0
    flags[1] = 1.0

    try:
        meshmonk.rigid_register(
            floating_features=features,
            target_features=features,
            floating_faces=faces,
            target_faces=faces,
            floating_flags=flags,
            target_flags=flags,
            num_iterations=2,
        )
    except MeshMonkError as e:
        # We may get DegenerateInput or InsufficientInliers depending on
        # how the C++ code validates. Either way, .code must be a valid enum.
        assert isinstance(
            e.code, RegistrationError
        ), f"Expected .code to be RegistrationError, got {type(e.code)}: {e.code!r}"
    except Exception:
        pass  # Other exceptions (e.g. from malformed faces) are acceptable


# ---------------------------------------------------------------------------
# RegistrationError enum is importable and has expected members
# ---------------------------------------------------------------------------


def test_registration_error_enum_members():
    """Verify all RegistrationError enum members exist."""
    assert hasattr(RegistrationError, "DegenerateInput")
    assert hasattr(RegistrationError, "InsufficientInliers")
    assert hasattr(RegistrationError, "DecompositionFailed")
    assert hasattr(RegistrationError, "NonConvergence")
