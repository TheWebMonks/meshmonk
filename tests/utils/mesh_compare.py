"""
Mesh comparison utilities for MeshMonk test harness.

Module constants (to be calibrated in v0.1):
    DEFAULT_TRANSFORM_ATOL: float — tolerance for rigid transform components
    DEFAULT_VERTEX_ATOL_MM: float — per-vertex position tolerance in millimetres
"""

import numpy as np
from scipy.spatial import cKDTree

# Placeholder tolerances — calibrate in v0.1 once run-to-run noise is measured
DEFAULT_TRANSFORM_ATOL = 1e-4
DEFAULT_VERTEX_ATOL_MM = 0.01


def _validate_vertex_array(arr: np.ndarray, name: str) -> np.ndarray:
    """Cast to float64 and validate shape (N, 3). Raises ValueError on bad input."""
    arr = np.asarray(arr, dtype=np.float64)
    if arr.shape[0] == 0:
        raise ValueError(f"empty input array: '{name}' has 0 rows")
    if arr.ndim != 2 or arr.shape[1] != 3:
        raise ValueError(
            f"invalid shape for '{name}': expected (N, 3), got {arr.shape}"
        )
    return arr


def rmse(a: np.ndarray, b: np.ndarray) -> float:
    """Root-mean-square per-vertex error between two (N, 3) vertex arrays.

    Both inputs are cast to float64 internally.
    Inputs must have the same shape (N, 3).
    Raises ValueError on empty or non-(N,3) input.
    """
    a = _validate_vertex_array(a, "a")
    b = _validate_vertex_array(b, "b")
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def max_vertex_distance(a: np.ndarray, b: np.ndarray) -> float:
    """L-infinity per-vertex distance between two (N, 3) vertex arrays.

    Both inputs are cast to float64 internally.
    Inputs must have the same shape (N, 3).
    Raises ValueError on empty or non-(N,3) input.
    """
    a = _validate_vertex_array(a, "a")
    b = _validate_vertex_array(b, "b")
    return float(np.max(np.sqrt(np.sum((a - b) ** 2, axis=1))))


def hausdorff_symmetric(a: np.ndarray, b: np.ndarray) -> float:
    """Symmetric Hausdorff distance between two (N, 3) and (M, 3) vertex arrays.

    Bidirectional: max(directed_hausdorff(a->b), directed_hausdorff(b->a)).
    Uses scipy.spatial.cKDTree for efficiency.
    Both inputs are cast to float64 internally.
    Raises ValueError on empty or non-(N,3) input.
    """
    a = _validate_vertex_array(a, "a")
    b = _validate_vertex_array(b, "b")
    tree_b = cKDTree(b)
    tree_a = cKDTree(a)
    dist_ab, _ = tree_b.query(a)
    dist_ba, _ = tree_a.query(b)
    return float(max(dist_ab.max(), dist_ba.max()))
