"""Benchmark: pyramid (multi-resolution nonrigid) registration at three mesh-size tiers.

Pyramid uses annealing so iteration count can vary run-to-run; the CI
regression tolerance for pyramid is therefore 30% (vs 20% for rigid/nonrigid).

Inputs are pre-aligned by a rigid pass (not timed). Pyramid runs 2 layers
in benchmarks (vs 3 in e2e tests) to keep wall-clock manageable.

Mesh tiers: 1K, 3K, 7K vertices (see benchmarks/README.md for rationale).
"""

from __future__ import annotations

import numpy as np
import pytest

import meshmonk
from benchmarks._harness import build_prealigned_inputs

MESH_SIZES = [1000, 3000, 7000]


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
    return build_prealigned_inputs(target_n)


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
