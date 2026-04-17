import tempfile
from pathlib import Path

import numpy as np
import pytest
from tests.utils.obj_to_npz import convert


MINIMAL_OBJ = """\
v 0.0 0.0 0.0
v 1.0 0.0 0.0
v 0.0 1.0 0.0
v 0.0 0.0 1.0
f 1 2 3
f 1 2 4
f 1 3 4
f 2 3 4
"""


def test_roundtrip_schema():
    with tempfile.TemporaryDirectory() as tmpdir:
        obj_path = Path(tmpdir) / "test.obj"
        npz_path = Path(tmpdir) / "test.npz"
        obj_path.write_text(MINIMAL_OBJ)

        convert(obj_path, npz_path)

        data = np.load(npz_path)
        assert "vertices" in data
        assert "faces" in data
        assert data["vertices"].shape == (4, 3)
        assert data["faces"].shape == (4, 3)
        assert data["vertices"].dtype == np.float64
        assert data["faces"].dtype == np.int32


def test_roundtrip_values():
    with tempfile.TemporaryDirectory() as tmpdir:
        obj_path = Path(tmpdir) / "test.obj"
        npz_path = Path(tmpdir) / "test.npz"
        obj_path.write_text(MINIMAL_OBJ)

        convert(obj_path, npz_path)

        data = np.load(npz_path)
        # Vertex (0,0,0) must be present
        assert np.any(np.all(data["vertices"] == [0.0, 0.0, 0.0], axis=1))
