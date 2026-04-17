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


def test_max_vertex_distance_known():
    """max_vertex_distance returns the L-inf per-vertex distance for a known case."""
    a = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    b = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    assert max_vertex_distance(a, b) == pytest.approx(1.0)


def test_hausdorff_known_value():
    """hausdorff_symmetric returns the correct value for a known geometric case."""
    a = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    b = np.array([[0.5, 0.0, 0.0]])
    # directed a->b: max distance from a's points to nearest in b: max(0.5, 0.5) = 0.5
    # directed b->a: distance from b's point to nearest in a: 0.5
    # symmetric: max(0.5, 0.5) = 0.5
    assert hausdorff_symmetric(a, b) == pytest.approx(0.5)


def test_rmse_empty_input_raises():
    """rmse raises ValueError for empty input arrays."""
    import pytest
    a = np.zeros((0, 3))
    b = np.zeros((0, 3))
    with pytest.raises(ValueError, match="empty input"):
        rmse(a, b)


def test_max_vertex_distance_empty_input_raises():
    """max_vertex_distance raises ValueError for empty input arrays."""
    import pytest
    a = np.zeros((0, 3))
    b = np.zeros((0, 3))
    with pytest.raises(ValueError, match="empty input"):
        max_vertex_distance(a, b)


def test_hausdorff_empty_input_raises():
    """hausdorff_symmetric raises ValueError for empty input arrays."""
    import pytest
    a = np.zeros((0, 3))
    b = np.zeros((0, 3))
    with pytest.raises(ValueError, match="empty input"):
        hausdorff_symmetric(a, b)


def test_rmse_wrong_ndim_raises():
    """rmse raises ValueError for non-(N,3) shaped input."""
    import pytest
    a = np.array([[0.0, 0.0], [1.0, 0.0]])  # (N, 2) — wrong
    b = np.array([[0.0, 0.0], [1.0, 0.0]])
    with pytest.raises(ValueError, match="shape"):
        rmse(a, b)
