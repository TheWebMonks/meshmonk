"""Shared harness utilities for MeshMonk benchmark suite.

Single canonical implementation of mesh-loading and downsampling helpers
used by bench_rigid.py, bench_nonrigid.py, bench_pyramid.py, and
test_harness_sizing.py.

Root-cause note (bead 03i.1)
-----------------------------
``meshmonk.downsample_mesh(ratio)`` interprets ``ratio`` as the **fraction of
vertices to REMOVE** (i.e. OpenMesh runs ``round(ratio * n_verts)`` edge
collapses), NOT the fraction to keep.

The original bug computed ``ratio = target_n / current_n``, which asked the
decimater to perform ``target_n`` collapses — removing almost all vertices
when the ratio was close to 1 (e.g. 7000/7160 ≈ 0.978 → removes 7000 verts,
leaving ~160).

The correct formula:

    ratio_to_remove = 1.0 - target_n / current_n

clamped to a maximum of 0.5 per step (never remove more than half the mesh
in a single pass) so we never overshoot on a large single decimation.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import meshmonk

DATA_DIR = Path(__file__).parent.parent / "data"
TEMPLATE_OBJ = DATA_DIR / "Template.obj"
DEMOFACE_OBJ = DATA_DIR / "demoFace.obj"


def load_full_mesh():
    """Load Template.obj and demoFace.obj as numpy arrays.

    Returns ``(float_v, float_f, target_v, target_f)`` — all float32/int32.
    Skips the calling test/benchmark if trimesh is not available.
    """
    trimesh = pytest.importorskip("trimesh")
    floating = trimesh.load(str(TEMPLATE_OBJ), process=False)
    target = trimesh.load(str(DEMOFACE_OBJ), process=False)
    float_v = np.asarray(floating.vertices, dtype=np.float32)
    float_f = np.asarray(floating.faces, dtype=np.int32)
    target_v = np.asarray(target.vertices, dtype=np.float32)
    target_f = np.asarray(target.faces, dtype=np.int32)
    return float_v, float_f, target_v, target_f


def downsample_to(feat, faces, flags, target_n: int):
    """Iteratively downsample a mesh to at most ``target_n`` vertices.

    ``meshmonk.downsample_mesh(ratio)`` interprets ``ratio`` as the fraction
    of vertices to REMOVE (number of edge collapses / current vertex count).
    The correct formula to reach ``target_n`` from ``current_n`` is:

        ratio_to_remove = 1.0 - target_n / current_n

    clamped to 0.5 maximum per step (never remove more than half the mesh
    in a single pass).

    The loop terminates when either:

    - current vertex count <= target_n, or
    - a decimation step makes no progress (boundary lock, degenerate mesh,
      or target already at minimum decimatable size).

    Parameters
    ----------
    feat:
        (N, 6) float32 feature matrix.
    faces:
        (M, 3) int32 face index array.
    flags:
        (N,) float32 per-vertex flags.
    target_n:
        Target vertex count (upper bound — actual result may be <= target_n).

    Returns
    -------
    (feat, faces, flags)
        Downsampled arrays with vertex count <= target_n.
    """
    current_n = feat.shape[0]
    while current_n > target_n:
        # Fraction of vertices to REMOVE — NOT the fraction to keep.
        ratio = 1.0 - target_n / current_n
        # Cap: never remove more than 50% per step to avoid overshoot.
        ratio = min(ratio, 0.5)
        feat, faces, flags, _ = meshmonk.downsample_mesh(feat, faces, flags, ratio)
        new_n = feat.shape[0]
        if new_n >= current_n:
            # Decimation made no progress (e.g. all boundary vertices locked).
            break
        current_n = new_n
    return feat, faces, flags
