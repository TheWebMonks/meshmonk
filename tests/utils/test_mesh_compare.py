import numpy as np
import pytest
from tests.utils.mesh_compare import (
    DEFAULT_TRANSFORM_ATOL,
    DEFAULT_VERTEX_ATOL_MM,
    hausdorff_symmetric,
    max_vertex_distance,
    rmse,
)


def test_rmse_identity():
    a = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    assert rmse(a, a) == pytest.approx(0.0)


def test_max_vertex_distance_identity():
    a = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    assert max_vertex_distance(a, a) == pytest.approx(0.0)


def test_hausdorff_identity():
    a = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    assert hausdorff_symmetric(a, a) == pytest.approx(0.0)


def test_rmse_known_translation():
    # Translate by (1, 0, 0): per-vertex distance is always 1.0, RMSE == 1.0
    N = 10
    a = np.zeros((N, 3))
    b = a.copy()
    b[:, 0] = 1.0
    assert rmse(a, b) == pytest.approx(1.0)


def test_hausdorff_symmetric_property():
    a = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    b = np.array([[0.5, 0.0, 0.0]])
    # hausdorff(a, b) == hausdorff(b, a)
    assert hausdorff_symmetric(a, b) == pytest.approx(hausdorff_symmetric(b, a))


def test_constants_exist():
    assert DEFAULT_TRANSFORM_ATOL == 1e-4
    assert DEFAULT_VERTEX_ATOL_MM == 0.01
