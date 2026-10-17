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


def test_meshmonk_error_has_code_property(empty_registration_kwargs):
    """MeshMonkError must have a .code attribute."""
    with pytest.raises(MeshMonkError) as exc_info:
        meshmonk.rigid_register(**empty_registration_kwargs)
    assert hasattr(
        exc_info.value, "code"
    ), "MeshMonkError instance must have a .code attribute"


def test_degenerate_input_code_is_registration_error_enum(empty_registration_kwargs):
    """MeshMonkError.code must return a RegistrationError enum member."""
    with pytest.raises(MeshMonkError) as exc_info:
        meshmonk.rigid_register(**empty_registration_kwargs)
    code = exc_info.value.code
    assert isinstance(
        code, RegistrationError
    ), f"Expected .code to be RegistrationError, got {type(code)}: {code!r}"


def test_degenerate_input_code_value(empty_registration_kwargs):
    """Empty input must produce DegenerateInput code."""
    with pytest.raises(MeshMonkError) as exc_info:
        meshmonk.rigid_register(**empty_registration_kwargs)
    assert exc_info.value.code == RegistrationError.DegenerateInput


# ---------------------------------------------------------------------------
# Programmatic error handling (the main use case)
# ---------------------------------------------------------------------------


def test_programmatic_error_dispatch(empty_registration_kwargs):
    """Users must be able to switch on .code instead of parsing strings."""
    try:
        meshmonk.rigid_register(**empty_registration_kwargs)
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


def test_except_runtime_error_catches_meshmonk_error(empty_registration_kwargs):
    """except RuntimeError must still catch MeshMonkError."""
    with pytest.raises(RuntimeError):
        meshmonk.rigid_register(**empty_registration_kwargs)


def test_except_meshmonk_error_catches_it(empty_registration_kwargs):
    """except MeshMonkError must catch the exception."""
    caught = False
    try:
        meshmonk.rigid_register(**empty_registration_kwargs)
    except MeshMonkError:
        caught = True
    assert caught


# ---------------------------------------------------------------------------
# Message string preserved
# ---------------------------------------------------------------------------


def test_error_message_still_contains_description(empty_registration_kwargs):
    """str(e) must still contain the human-readable error description."""
    with pytest.raises(MeshMonkError) as exc_info:
        meshmonk.rigid_register(**empty_registration_kwargs)
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
    features = np.random.RandomState(0).randn(n, 6).astype(np.float32)
    faces = np.zeros((1, 3), dtype=np.int32)  # minimal faces
    # Set all flags to near-zero except 2 (fewer than the 4 minimum)
    flags = np.zeros(n, dtype=np.float32)
    flags[0] = 1.0
    flags[1] = 1.0

    with pytest.raises(MeshMonkError) as exc_info:
        meshmonk.rigid_register(
            floating_features=features,
            target_features=features,
            floating_faces=faces,
            target_faces=faces,
            floating_flags=flags,
            target_flags=flags,
            num_iterations=0,  # straight to the final inlier check (iterating can fail earlier in the eigensolver)
        )
    assert exc_info.value.code == RegistrationError.InsufficientInliers


# ---------------------------------------------------------------------------
# RegistrationError enum is importable and has expected members
# ---------------------------------------------------------------------------


def test_registration_error_enum_members():
    """Verify all RegistrationError enum members exist."""
    assert hasattr(RegistrationError, "DegenerateInput")
    assert hasattr(RegistrationError, "InsufficientInliers")
    assert hasattr(RegistrationError, "DecompositionFailed")
    assert hasattr(RegistrationError, "NonConvergence")
