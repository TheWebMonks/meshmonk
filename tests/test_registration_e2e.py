"""End-to-end registration tests using real Template.obj + demoFace.obj.

Mirrors the demo scripts under demo/ (rigid, nonrigid, pyramid, full pipeline)
and asserts numeric invariants so refactors are guarded against silent
regressions. Total runtime is ~6 s — unmarked so they run by default.

Skipped if the data files are missing or trimesh is not installed.
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
    """Load Template.obj (floating) and demoFace.obj (target) once per module."""
    trimesh = pytest.importorskip("trimesh")
    if not TEMPLATE_OBJ.exists() or not DEMOFACE_OBJ.exists():
        pytest.skip(f"OBJ data missing in {DATA_DIR}")

    floating = trimesh.load(str(TEMPLATE_OBJ), process=False)
    target = trimesh.load(str(DEMOFACE_OBJ), process=False)

    return {
        "floating_v": np.asarray(floating.vertices, dtype=np.float32),
        "floating_f": np.asarray(floating.faces, dtype=np.int32),
        "target_v": np.asarray(target.vertices, dtype=np.float32),
        "target_f": np.asarray(target.faces, dtype=np.int32),
    }


@pytest.fixture(scope="module")
def features(meshes):
    feat_float = meshmonk.features_from_vertices(
        meshes["floating_v"], meshes["floating_f"]
    )
    feat_target = meshmonk.features_from_vertices(
        meshes["target_v"], meshes["target_f"]
    )
    return {"float": feat_float, "target": feat_target}


@pytest.fixture(scope="module")
def rigid_result(meshes, features):
    """Run rigid once; reused by nonrigid / pyramid tests as their input."""
    meshmonk.set_log_level("warning")
    return meshmonk.rigid_register(
        floating_features=features["float"],
        target_features=features["target"],
        floating_faces=meshes["floating_f"],
        target_faces=meshes["target_f"],
        num_iterations=40,
        use_scaling=True,
    )


# ---------------------------------------------------------------------------
# Rigid
# ---------------------------------------------------------------------------


def test_rigid_produces_valid_alignment(meshes, rigid_result):
    r = rigid_result
    n = meshes["floating_v"].shape[0]

    assert r.aligned_vertices.shape == (n, 3)
    assert r.aligned_features.shape == (n, 6)
    assert np.isfinite(r.aligned_vertices).all()
    assert r.iterations_run > 0

    # Rigid moves the template a lot (~60 units) onto the target.
    move = np.linalg.norm(r.aligned_vertices - meshes["floating_v"], axis=1)
    assert move.mean() > 10.0, f"rigid mean-move unexpectedly small: {move.mean():.2f}"

    # Aligned template centroid should sit inside the target's bounding box
    # (template is just a face; target is a full head scan, so centroids
    # don't coincide — but the face must land on the head).
    aligned_c = r.aligned_vertices.mean(axis=0)
    tlo = meshes["target_v"].min(axis=0)
    thi = meshes["target_v"].max(axis=0)
    assert np.all(aligned_c >= tlo) and np.all(
        aligned_c <= thi
    ), f"rigid output centroid {aligned_c} outside target bbox [{tlo}, {thi}]"


# ---------------------------------------------------------------------------
# Nonrigid (viscoelastic, single-resolution) — on rigid output
# ---------------------------------------------------------------------------


def test_nonrigid_refines_rigid_output(meshes, features, rigid_result):
    r = meshmonk.nonrigid_register(
        floating_features=rigid_result.aligned_features,
        target_features=features["target"],
        floating_faces=meshes["floating_f"],
        target_faces=meshes["target_f"],
        num_iterations=40,
    )

    n = meshes["floating_v"].shape[0]
    assert r.aligned_vertices.shape == (n, 3)
    assert np.isfinite(r.aligned_vertices).all()
    assert r.final_inlier_weights.shape == (n,)
    assert r.final_inlier_weights.mean() > 0.5

    # Nonrigid should make small-but-nonzero corrections on top of rigid.
    step = np.linalg.norm(r.aligned_vertices - rigid_result.aligned_vertices, axis=1)
    assert (
        0.1 < step.mean() < 10.0
    ), f"nonrigid step-mean out of range: {step.mean():.2f}"
    assert step.max() < 50.0


# ---------------------------------------------------------------------------
# Pyramid (viscoelastic, multi-resolution) — on rigid output
# ---------------------------------------------------------------------------


def test_pyramid_refines_rigid_output(meshes, features, rigid_result):
    r = meshmonk.pyramid_register(
        floating_features=rigid_result.aligned_features,
        target_features=features["target"],
        floating_faces=meshes["floating_f"],
        target_faces=meshes["target_f"],
        num_iterations=10,
        num_pyramid_layers=3,
    )

    n = meshes["floating_v"].shape[0]
    assert r.aligned_vertices.shape == (n, 3)
    assert np.isfinite(r.aligned_vertices).all()
    assert r.final_inlier_weights.mean() > 0.5
    assert len(r.per_layer_iterations) == 3
    assert all(i > 0 for i in r.per_layer_iterations)

    step = np.linalg.norm(r.aligned_vertices - rigid_result.aligned_vertices, axis=1)
    assert (
        0.1 < step.mean() < 10.0
    ), f"pyramid step-mean out of range: {step.mean():.2f}"
    assert step.max() < 50.0


# ---------------------------------------------------------------------------
# Full pipeline: rigid, then nonrigid vs pyramid as alternatives
# ---------------------------------------------------------------------------


def test_e2e_pipeline_consistency(meshes, features, rigid_result):
    """Sanity check that nonrigid and pyramid, as alternatives on rigid output,
    both produce similarly-shaped results (not identical, but in the same ballpark)."""
    nr = meshmonk.nonrigid_register(
        floating_features=rigid_result.aligned_features,
        target_features=features["target"],
        floating_faces=meshes["floating_f"],
        target_faces=meshes["target_f"],
        num_iterations=40,
    )
    py = meshmonk.pyramid_register(
        floating_features=rigid_result.aligned_features,
        target_features=features["target"],
        floating_faces=meshes["floating_f"],
        target_faces=meshes["target_f"],
        num_iterations=10,
        num_pyramid_layers=3,
    )

    # The two methods should agree to within a few units per vertex on average.
    diff = np.linalg.norm(nr.aligned_vertices - py.aligned_vertices, axis=1)
    assert (
        diff.mean() < 5.0
    ), f"nonrigid vs pyramid mean-divergence too high: {diff.mean():.2f}"
