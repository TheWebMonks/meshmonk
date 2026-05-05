"""Memory-layout regression guard for InlierDetector, NonrigidRegistration, ViscoElasticTransformer.

Asserts that `nonrigid_register` output matches a committed reference within a
cross-platform tolerance, to catch large accidental changes to operation order,
memory layout, or numerical behavior in those classes.

Bitwise equality (np.array_equal) is not used: FP non-associativity, FMA, and
math-library differences cause genuine cross-platform drift even on the same
inputs (observed up to ~2.6e-2 between Linux x86_64 and macOS arm64). The
tolerance is loose enough to ride that variance but tight enough to catch an
order-of-magnitude algorithm regression. For strict same-process equivalence,
see `test_back_to_back_calls_same_size`.
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


def test_close_to_reference(registration_result):
    """Cross-platform regression guard against committed reference output.

    If this fails, operations on InlierDetector, NonrigidRegistration, or
    ViscoElasticTransformer changed by more than the cross-platform FP drift
    budget. Investigate whether the change is intentional.

    Tolerance chosen empirically from observed cross-platform drift on the
    same source: ~3e-4 macOS-arm64, ~1e-2 Windows-x86_64, ~3e-2 Linux-x86_64
    (all relative to a macOS-captured reference).
    """
    reference = np.load(str(REFERENCE_NPY))
    output = registration_result

    assert output.shape == reference.shape, (
        f"Shape mismatch: {output.shape} vs {reference.shape}"
    )
    np.testing.assert_allclose(
        output,
        reference,
        rtol=5e-2,
        atol=5e-2,
        err_msg=(
            "nonrigid_register output drifted from reference beyond the "
            "cross-platform FP budget — likely a real algorithm regression "
            "in InlierDetector / NonrigidRegistration / ViscoElasticTransformer."
        ),
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
