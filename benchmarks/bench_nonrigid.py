"""Benchmark: nonrigid (viscoelastic) registration at three mesh-size tiers.

Inputs are pre-aligned by a rigid pass so nonrigid has realistic starting
conditions. The rigid pre-alignment is done once per mesh-size tier in the
module-scoped fixture and is NOT included in the benchmark timer.

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
def nonrigid_inputs(request):
    """Pre-built inputs (post-rigid) for nonrigid registration."""
    target_n = request.param
    return build_prealigned_inputs(target_n)


# ---------------------------------------------------------------------------
# Benchmarks
# ---------------------------------------------------------------------------


def test_bench_nonrigid(benchmark, nonrigid_inputs):
    """Wall-clock time for nonrigid_register() at each mesh-size tier."""
    aligned_feat, float_f, target_feat, target_f = nonrigid_inputs

    meshmonk.set_log_level("silent")

    result = benchmark.pedantic(
        meshmonk.nonrigid_register,
        kwargs={
            "floating_features": aligned_feat,
            "target_features": target_feat,
            "floating_faces": float_f,
            "target_faces": target_f,
            "num_iterations": 20,
        },
        rounds=5,
        warmup_rounds=1,
    )

    assert result.aligned_vertices.shape == (aligned_feat.shape[0], 3)
    assert np.isfinite(result.aligned_vertices).all()
