import pytest
import meshmonk
import numpy as np


def test_meshmonk_error_is_runtime_error_subclass():
    assert issubclass(meshmonk.MeshMonkError, RuntimeError)


def test_degenerate_input_raises_meshmonk_error(empty_registration_kwargs):
    with pytest.raises(meshmonk.MeshMonkError):
        meshmonk.rigid_register(**empty_registration_kwargs)


def test_meshmonk_error_message_identifies_error_type(empty_registration_kwargs):
    with pytest.raises(meshmonk.MeshMonkError) as exc_info:
        meshmonk.rigid_register(**empty_registration_kwargs)
    assert "DegenerateInput" in str(exc_info.value)


def test_meshmonk_error_is_not_same_as_runtime_error():
    assert meshmonk.MeshMonkError is not RuntimeError
    assert issubclass(meshmonk.MeshMonkError, RuntimeError)
