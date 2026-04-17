import tempfile
from pathlib import Path

import numpy as np
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
        # All face indices must be present (0-based after OBJ 1-based conversion)
        expected_faces = sorted([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
        assert sorted(data["faces"].tolist()) == expected_faces


MINIMAL_OBJ_WITH_NORMALS = """\
v 0.0 0.0 0.0
v 1.0 0.0 0.0
v 0.0 1.0 0.0
v 0.0 0.0 1.0
vn 0.0 0.0 1.0
vn 1.0 0.0 0.0
vn 0.0 1.0 0.0
vn 0.0 0.0 1.0
f 1//1 2//2 3//3
f 1//1 2//2 4//4
f 1//1 3//3 4//4
f 2//2 3//3 4//4
"""


def test_roundtrip_normals():
    """OBJ vertex normals (obj:vn key in meshio) survive the round-trip."""
    with tempfile.TemporaryDirectory() as tmpdir:
        obj_path = Path(tmpdir) / "test_normals.obj"
        npz_path = Path(tmpdir) / "test_normals.npz"
        obj_path.write_text(MINIMAL_OBJ_WITH_NORMALS)

        convert(obj_path, npz_path)

        data = np.load(npz_path)
        assert "normals" in data, "normals key missing from NPZ"
        assert data["normals"].shape == (4, 3)
        assert data["normals"].dtype == np.float64


def test_non_triangle_cells_raise():
    """convert() raises ValueError when non-triangle cells are present alongside triangles."""
    quad_obj = """\
v 0.0 0.0 0.0
v 1.0 0.0 0.0
v 1.0 1.0 0.0
v 0.0 1.0 0.0
v 0.0 0.0 1.0
f 1 2 3
f 1 2 3 4
"""
    with tempfile.TemporaryDirectory() as tmpdir:
        obj_path = Path(tmpdir) / "mixed.obj"
        npz_path = Path(tmpdir) / "mixed.npz"
        obj_path.write_text(quad_obj)

        import pytest
        with pytest.raises(ValueError, match="non-triangle cells"):
            convert(obj_path, npz_path)


def test_arrays_are_c_contiguous():
    """Vertices and faces must be C-contiguous for downstream Eigen::Map consumers."""
    with tempfile.TemporaryDirectory() as tmpdir:
        obj_path = Path(tmpdir) / "test.obj"
        npz_path = Path(tmpdir) / "test.npz"
        obj_path.write_text(MINIMAL_OBJ)

        convert(obj_path, npz_path)

        data = np.load(npz_path)
        assert data["vertices"].flags["C_CONTIGUOUS"]
        assert data["faces"].flags["C_CONTIGUOUS"]
