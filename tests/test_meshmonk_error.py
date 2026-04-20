import pytest
import meshmonk
import numpy as np


def test_meshmonk_error_is_runtime_error_subclass():
    assert issubclass(meshmonk.MeshMonkError, RuntimeError)


def test_degenerate_input_raises_meshmonk_error():
    empty = np.zeros((0, 6), dtype=np.float32)
    faces = np.zeros((0, 3), dtype=np.int32)
    flags = np.zeros((0,), dtype=np.float32)
    with pytest.raises(meshmonk.MeshMonkError):
        meshmonk.rigid_register(
            floating_features=empty,
            target_features=empty,
            floating_faces=faces,
            target_faces=faces,
            floating_flags=flags,
            target_flags=flags,
        )


def test_meshmonk_error_message_identifies_error_type():
    empty = np.zeros((0, 6), dtype=np.float32)
    faces = np.zeros((0, 3), dtype=np.int32)
    flags = np.zeros((0,), dtype=np.float32)
    with pytest.raises(meshmonk.MeshMonkError) as exc_info:
        meshmonk.rigid_register(
            floating_features=empty,
            target_features=empty,
            floating_faces=faces,
            target_faces=faces,
            floating_flags=flags,
            target_flags=flags,
        )
    assert "DegenerateInput" in str(exc_info.value)


def test_meshmonk_error_is_not_same_as_runtime_error():
    assert meshmonk.MeshMonkError is not RuntimeError
    assert issubclass(meshmonk.MeshMonkError, RuntimeError)
