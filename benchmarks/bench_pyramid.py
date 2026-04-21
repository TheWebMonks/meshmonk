"""Benchmark: pyramid (multi-resolution nonrigid) registration at three mesh-size tiers.

Pyramid uses annealing so iteration count can vary run-to-run; the CI
regression tolerance for pyramid is therefore 30% (vs 20% for rigid/nonrigid).

Inputs are pre-aligned by a rigid pass (not timed). Pyramid runs 2 layers
in benchmarks (vs 3 in e2e tests) to keep wall-clock manageable.

Mesh tiers: 1K, 3K, 7K vertices (see benchmarks/README.md for rationale).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import meshmonk

DATA_DIR = Path(__file__).parent.parent / "data"
TEMPLATE_OBJ = DATA_DIR / "Template.obj"
DEMOFACE_OBJ = DATA_DIR / "demoFace.obj"

MESH_SIZES = [1000, 3000, 7000]
RNG_SEED = 42


def _load_full_mesh():
    trimesh = pytest.importorskip("trimesh")
    floating = trimesh.load(str(TEMPLATE_OBJ), process=False)
    target = trimesh.load(str(DEMOFACE_OBJ), process=False)
    float_v = np.asarray(floating.vertices, dtype=np.float32)
    float_f = np.asarray(floating.faces, dtype=np.int32)
    target_v = np.asarray(target.vertices, dtype=np.float32)
    target_f = np.asarray(target.faces, dtype=np.int32)
    return float_v, float_f, target_v, target_f


def _downsample_to(feat, faces, flags, target_n: int):
    current_n = feat.shape[0]
    while current_n > target_n:
        ratio = target_n / current_n
        ratio = max(ratio, 0.5)
        feat, faces, flags, _ = meshmonk.downsample_mesh(feat, faces, flags, ratio)
        new_n = feat.shape[0]
        if new_n >= current_n:
            break
        current_n = new_n
    return feat, faces, flags


def _build_pyramid_inputs(target_n: int):
    """Build inputs for pyramid: floating pre-aligned via rigid."""
    float_v, float_f, target_v, target_f = _load_full_mesh()

    full_n = float_v.shape[0]
    float_feat = meshmonk.features_from_vertices(float_v, float_f)
    float_flags = np.ones(full_n, dtype=np.float32)
    target_feat = meshmonk.features_from_vertices(target_v, target_f)

    if target_n < full_n:
        float_feat, float_f, float_flags = _downsample_to(
            float_feat, float_f, float_flags, target_n
        )

    meshmonk.set_log_level("silent")
    rigid_result = meshmonk.rigid_register(
        floating_features=float_feat,
        target_features=target_feat,
        floating_faces=float_f,
        target_faces=target_f,
        num_iterations=20,
        use_scaling=True,
    )
    aligned_feat = rigid_result.aligned_features

    return aligned_feat, float_f, target_feat, target_f


# ---------------------------------------------------------------------------
# Parametrised fixture
# ---------------------------------------------------------------------------


@pytest.fixture(
    params=MESH_SIZES,
    ids=[f"{n // 1000}K" for n in MESH_SIZES],
    scope="module",
)
def pyramid_inputs(request):
    """Pre-built inputs (post-rigid) for pyramid registration."""
    target_n = request.param
    return _build_pyramid_inputs(target_n)


# ---------------------------------------------------------------------------
# Benchmarks
# ---------------------------------------------------------------------------


def test_bench_pyramid(benchmark, pyramid_inputs):
    """Wall-clock time for pyramid_register() at each mesh-size tier."""
    aligned_feat, float_f, target_feat, target_f = pyramid_inputs

    meshmonk.set_log_level("silent")

    result = benchmark.pedantic(
        meshmonk.pyramid_register,
        kwargs={
            "floating_features": aligned_feat,
            "target_features": target_feat,
            "floating_faces": float_f,
            "target_faces": target_f,
            "num_iterations": 10,
            "num_pyramid_layers": 2,
        },
        rounds=5,
        warmup_rounds=1,
    )

    assert result.aligned_vertices.shape == (aligned_feat.shape[0], 3)
    assert np.isfinite(result.aligned_vertices).all()
    assert len(result.per_layer_iterations) == 2
