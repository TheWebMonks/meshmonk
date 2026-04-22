"""Memory-layout regression guard for InlierDetector, NonrigidRegistration, ViscoElasticTransformer.

Asserts bitwise equality of `nonrigid_register` output against a committed reference to catch
accidental changes to operation order, memory layout, or numerical behavior in those classes.

Use numpy.array_equal (bitwise), NOT numpy.testing.assert_allclose.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import meshmonk

DATA_DIR = Path(__file__).parent.parent / "data"
TEMPLATE_OBJ = DATA_DIR / "Template.obj"
DEMOFACE_OBJ = DATA_DIR / "demoFace.obj"
REFERENCE_NPY = Path(__file__).parent / "golden" / "memory_layout_reference.npy"


@pytest.fixture(scope="module")
def registration_result():
    """Run the same nonrigid registration used to capture the reference."""
    trimesh = pytest.importorskip("trimesh")
    if not TEMPLATE_OBJ.exists() or not DEMOFACE_OBJ.exists():
        pytest.skip(f"OBJ data missing in {DATA_DIR}")
    if not REFERENCE_NPY.exists():
        pytest.skip(f"Reference file missing: {REFERENCE_NPY}")

    floating = trimesh.load(str(TEMPLATE_OBJ), process=False)
    target = trimesh.load(str(DEMOFACE_OBJ), process=False)

    float_v = np.asarray(floating.vertices, dtype=np.float32)
    float_f = np.asarray(floating.faces, dtype=np.int32)
    target_v = np.asarray(target.vertices, dtype=np.float32)
    target_f = np.asarray(target.faces, dtype=np.int32)

    float_feat = meshmonk.features_from_vertices(float_v, float_f)
    target_feat = meshmonk.features_from_vertices(target_v, target_f)

    meshmonk.set_log_level("silent")

    result = meshmonk.nonrigid_register(
        floating_features=float_feat,
        target_features=target_feat,
        floating_faces=float_f,
        target_faces=target_f,
        num_iterations=5,
    )
    return result.aligned_vertices.copy()


def test_bitwise_equality_with_reference(registration_result):
    """Bitwise regression guard against committed reference output.

    If this fails, operations on InlierDetector, NonrigidRegistration, or
    ViscoElasticTransformer changed. Investigate whether the change is intentional.
    Use numpy.array_equal (bitwise), not assert_allclose.
    """
    reference = np.load(str(REFERENCE_NPY))
    output = registration_result

    assert output.shape == reference.shape, (
        f"Shape mismatch: {output.shape} vs {reference.shape}"
    )
    assert np.array_equal(output, reference), (
        "Bitwise equality FAILED — operations on InlierDetector / NonrigidRegistration / "
        "ViscoElasticTransformer changed computed values. "
        f"Max diff: {np.abs(output - reference).max():.2e}, "
        f"Mismatched elements: {(output != reference).sum()} / {reference.size}"
    )


def test_back_to_back_calls_same_size(registration_result):
    """Two consecutive nonrigid_register calls with same-size mesh give identical results.

    Guards against state leakage between calls to `nonrigid_register` — a broad correctness
    invariant for the Registration classes (InlierDetector, NonrigidRegistration,
    ViscoElasticTransformer).
    """
    trimesh = pytest.importorskip("trimesh")
    if not TEMPLATE_OBJ.exists() or not DEMOFACE_OBJ.exists():
        pytest.skip(f"OBJ data missing in {DATA_DIR}")

    floating = trimesh.load(str(TEMPLATE_OBJ), process=False)
    target = trimesh.load(str(DEMOFACE_OBJ), process=False)

    float_v = np.asarray(floating.vertices, dtype=np.float32)
    float_f = np.asarray(floating.faces, dtype=np.int32)
    target_v = np.asarray(target.vertices, dtype=np.float32)
    target_f = np.asarray(target.faces, dtype=np.int32)

    float_feat = meshmonk.features_from_vertices(float_v, float_f)
    target_feat = meshmonk.features_from_vertices(target_v, target_f)

    meshmonk.set_log_level("silent")

    r1 = meshmonk.nonrigid_register(
        floating_features=float_feat,
        target_features=target_feat,
        floating_faces=float_f,
        target_faces=target_f,
        num_iterations=5,
    )
    # Re-compute features (float_feat is mutated in-place by registration)
    float_feat2 = meshmonk.features_from_vertices(float_v, float_f)
    r2 = meshmonk.nonrigid_register(
        floating_features=float_feat2,
        target_features=target_feat,
        floating_faces=float_f,
        target_faces=target_f,
        num_iterations=5,
    )

    assert np.array_equal(r1.aligned_vertices, r2.aligned_vertices), (
        "Back-to-back calls gave different results — state leaked between calls to nonrigid_register."
    )
