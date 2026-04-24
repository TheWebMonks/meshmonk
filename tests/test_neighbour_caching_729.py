"""Correctness guard for the NeighbourFinder warm-start cache (bead 729).

The cache pre-seeds nanoflann's KNNResultSet worstDist with the previous
iteration's k-th distance at current positions, compensated by (1+eps) to
cancel the approximate-search sibling-descent multiplier. The acceptance
criterion for bead 729 requires that cached and uncached code paths
produce identical registration output — this test enforces that.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import meshmonk

DATA_DIR = Path(__file__).parent.parent / "data"
TEMPLATE_OBJ = DATA_DIR / "Template.obj"
DEMOFACE_OBJ = DATA_DIR / "demoFace.obj"


@pytest.fixture(scope="module")
def meshes():
    trimesh = pytest.importorskip("trimesh")
    if not TEMPLATE_OBJ.exists() or not DEMOFACE_OBJ.exists():
        pytest.skip(f"OBJ data missing in {DATA_DIR}")
    floating = trimesh.load(str(TEMPLATE_OBJ), process=False)
    target = trimesh.load(str(DEMOFACE_OBJ), process=False)
    return (
        np.asarray(floating.vertices, dtype=np.float32),
        np.asarray(floating.faces, dtype=np.int32),
        np.asarray(target.vertices, dtype=np.float32),
        np.asarray(target.faces, dtype=np.int32),
    )


def _register(float_v, float_f, target_v, target_f, iterations):
    float_feat = meshmonk.features_from_vertices(float_v, float_f)
    target_feat = meshmonk.features_from_vertices(target_v, target_f)
    return meshmonk.nonrigid_register(
        floating_features=float_feat,
        target_features=target_feat,
        floating_faces=float_f,
        target_faces=target_f,
        num_iterations=iterations,
    ).aligned_vertices.copy()


def test_cached_matches_uncached_nonrigid(meshes):
    """Registration output is bit-identical whether caching is on or off."""
    if not hasattr(meshmonk, "_set_neighbour_caching"):
        pytest.skip("build lacks _set_neighbour_caching (pre-729 build)")

    meshmonk.set_log_level("silent")
    float_v, float_f, target_v, target_f = meshes

    original = meshmonk._get_neighbour_caching()
    try:
        meshmonk._set_neighbour_caching(False)
        uncached = _register(float_v, float_f, target_v, target_f, iterations=5)

        meshmonk._set_neighbour_caching(True)
        cached = _register(float_v, float_f, target_v, target_f, iterations=5)
    finally:
        meshmonk._set_neighbour_caching(original)

    assert np.array_equal(cached, uncached), (
        "Cached and uncached paths diverged — max diff "
        f"{np.abs(cached - uncached).max():.3e}, "
        f"mismatched {(cached != uncached).sum()} / {cached.size}"
    )


def test_caching_default_is_on():
    if not hasattr(meshmonk, "_get_neighbour_caching"):
        pytest.skip("build lacks _get_neighbour_caching (pre-729 build)")
    assert meshmonk._get_neighbour_caching() is True
