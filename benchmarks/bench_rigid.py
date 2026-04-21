"""Benchmark: rigid registration at three mesh-size tiers.

Mesh tiers: 1K, 3K, 7K vertices.
  - 7K: full Template.obj (7,160 vertices) — used as-is.
  - 3K, 1K: obtained by iterative downsampling from the full mesh using
    meshmonk.downsample_mesh with a deterministic ratio.

All inputs are synthesised from data/Template.obj with a fixed random seed
perturbation so benchmarks are fully reproducible.

See benchmarks/README.md for tier-choice rationale.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import meshmonk

DATA_DIR = Path(__file__).parent.parent / "data"
TEMPLATE_OBJ = DATA_DIR / "Template.obj"
DEMOFACE_OBJ = DATA_DIR / "demoFace.obj"

# Mesh-size tiers (vertex counts). The top tier matches Template.obj (7,160 v).
MESH_SIZES = [1000, 3000, 7000]

# Fixed random state for reproducible perturbation
RNG_SEED = 42


def _load_full_mesh():
    """Load Template.obj and demoFace.obj as numpy arrays."""
    trimesh = pytest.importorskip("trimesh")
    floating = trimesh.load(str(TEMPLATE_OBJ), process=False)
    target = trimesh.load(str(DEMOFACE_OBJ), process=False)
    float_v = np.asarray(floating.vertices, dtype=np.float32)
    float_f = np.asarray(floating.faces, dtype=np.int32)
    target_v = np.asarray(target.vertices, dtype=np.float32)
    target_f = np.asarray(target.faces, dtype=np.int32)
    return float_v, float_f, target_v, target_f


def _downsample_to(feat, faces, flags, target_n: int):
    """Iteratively downsample until vertex count <= target_n."""
    current_n = feat.shape[0]
    while current_n > target_n:
        ratio = target_n / current_n
        # Clamp ratio to avoid over-aggressive single steps
        ratio = max(ratio, 0.5)
        feat, faces, flags, _ = meshmonk.downsample_mesh(feat, faces, flags, ratio)
        new_n = feat.shape[0]
        if new_n >= current_n:
            break
        current_n = new_n
    return feat, faces, flags


def _build_inputs(target_n: int):
    """Build (floating_feat, floating_faces, target_feat, target_faces) at target_n vertices."""
    float_v, float_f, target_v, target_f = _load_full_mesh()

    full_n = float_v.shape[0]
    float_feat = meshmonk.features_from_vertices(float_v, float_f)
    float_flags = np.ones(full_n, dtype=np.float32)

    target_feat = meshmonk.features_from_vertices(target_v, target_f)

    if target_n < full_n:
        float_feat, float_f, float_flags = _downsample_to(
            float_feat, float_f, float_flags, target_n
        )

    # Apply a small known perturbation to floating so rigid has something to do.
    # Fixed seed for reproducibility.
    rng = np.random.default_rng(RNG_SEED)
    perturb = rng.uniform(-2.0, 2.0, size=(float_feat.shape[0], 3)).astype(np.float32)
    float_feat = float_feat.copy()
    float_feat[:, :3] += perturb

    return float_feat, float_f, target_feat, target_f


# ---------------------------------------------------------------------------
# Parametrised fixture over mesh sizes
# ---------------------------------------------------------------------------


@pytest.fixture(
    params=MESH_SIZES,
    ids=[f"{n // 1000}K" for n in MESH_SIZES],
    scope="module",
)
def rigid_inputs(request):
    """Pre-built inputs for rigid registration at each mesh-size tier."""
    target_n = request.param
    return _build_inputs(target_n)


# ---------------------------------------------------------------------------
# Benchmarks
# ---------------------------------------------------------------------------


def test_bench_rigid(benchmark, rigid_inputs):
    """Wall-clock time for rigid_register() at each mesh-size tier."""
    float_feat, float_f, target_feat, target_f = rigid_inputs

    meshmonk.set_log_level("silent")

    result = benchmark.pedantic(
        meshmonk.rigid_register,
        kwargs={
            "floating_features": float_feat,
            "target_features": target_feat,
            "floating_faces": float_f,
            "target_faces": target_f,
            "num_iterations": 20,
            "use_scaling": True,
        },
        rounds=5,
        warmup_rounds=1,
    )

    # Sanity: result must be a valid alignment
    assert result.aligned_vertices.shape == (float_feat.shape[0], 3)
    assert np.isfinite(result.aligned_vertices).all()
