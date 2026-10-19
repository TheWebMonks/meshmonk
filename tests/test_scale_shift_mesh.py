"""
meshmonk.scale_shift_mesh — transfer registered features from a low-res pyramid
layer to a higher-res one, as pyramid_register does between layers
(PyramidNonrigidRegistration.cpp, ScaleShifter).

Matching nodes (original index present in both layers) take the low-res
features verbatim. New nodes start from their original high-res features and
move by the inverse-distance-weighted deformation of their matching
neighbours. Under a uniform translation every neighbour deformation is the
same vector, so each new node must land exactly on its original position plus
that translation: the hand-computed expectation below.
"""

import numpy as np
import pytest

import meshmonk

trimesh = pytest.importorskip("trimesh")

TRANSLATION = np.array([1.5, -2.0, 0.75], dtype=np.float32)


@pytest.fixture(scope="module")
def pyramid_pair():
    m = trimesh.creation.icosphere(subdivisions=3)
    high_features = np.hstack([m.vertices, m.vertex_normals]).astype(np.float32)
    faces = np.asarray(m.faces, dtype=np.int32)
    flags = np.ones(len(high_features), dtype=np.float32)
    low_features, _, _, low_indices = meshmonk.downsample_mesh(
        high_features, faces, flags, downsample_ratio=0.5
    )
    high_indices = np.arange(len(high_features), dtype=np.int32)
    return high_features, high_indices, np.asarray(low_features), low_indices


def _shift(pair, previous_features):
    high_features, high_indices, _, low_indices = pair
    return meshmonk.scale_shift_mesh(
        previous_features, low_indices, high_features, high_indices
    )


def test_matching_nodes_take_low_res_features(pyramid_pair):
    high_features, _, low_features, low_indices = pyramid_pair
    registered = low_features.copy()
    registered[:, :3] += TRANSLATION

    out = _shift(pyramid_pair, registered)

    assert out.shape == high_features.shape
    np.testing.assert_array_equal(out[low_indices], registered)


def test_new_nodes_follow_uniform_translation(pyramid_pair):
    high_features, _, low_features, low_indices = pyramid_pair
    registered = low_features.copy()
    registered[:, :3] += TRANSLATION

    out = _shift(pyramid_pair, registered)

    new = np.setdiff1d(np.arange(len(high_features)), low_indices)
    assert len(new) > 0
    assert np.isfinite(out).all()
    np.testing.assert_allclose(
        out[new, :3], high_features[new, :3] + TRANSLATION, atol=1e-5
    )
    # ScaleShifter moves positions only; new nodes keep their original normals.
    np.testing.assert_array_equal(out[new, 3:], high_features[new, 3:])


def test_rejects_mismatched_shapes(pyramid_pair):
    high_features, high_indices, low_features, low_indices = pyramid_pair
    with pytest.raises(ValueError, match="new_features"):
        meshmonk.scale_shift_mesh(
            low_features, low_indices, high_features[:-1], high_indices
        )
    with pytest.raises(ValueError, match="previous_features"):
        meshmonk.scale_shift_mesh(
            low_features[:-1], low_indices, high_features, high_indices
        )
