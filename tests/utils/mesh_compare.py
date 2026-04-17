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


def rmse(a: np.ndarray, b: np.ndarray) -> float:
    """Root-mean-square per-vertex error between two (N, 3) vertex arrays.

    Both inputs are cast to float64 internally.
    Inputs must have the same shape (N, 3).
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def max_vertex_distance(a: np.ndarray, b: np.ndarray) -> float:
    """L-infinity per-vertex distance between two (N, 3) vertex arrays.

    Both inputs are cast to float64 internally.
    Inputs must have the same shape (N, 3).
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    return float(np.max(np.sqrt(np.sum((a - b) ** 2, axis=1))))


def hausdorff_symmetric(a: np.ndarray, b: np.ndarray) -> float:
    """Symmetric Hausdorff distance between two (N, 3) and (M, 3) vertex arrays.

    Bidirectional: max(directed_hausdorff(a->b), directed_hausdorff(b->a)).
    Uses scipy.spatial.cKDTree for efficiency.
    Both inputs are cast to float64 internally.
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    tree_b = cKDTree(b)
    tree_a = cKDTree(a)
    dist_ab, _ = tree_b.query(a)
    dist_ba, _ = tree_a.query(b)
    return float(max(dist_ab.max(), dist_ba.max()))
