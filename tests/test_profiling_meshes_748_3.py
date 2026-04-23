"""Tests for profiling mesh generation script (bead meshmonk-modernization-748.3).

Verifies:
1. The generate_profiling_meshes.py script produces 6 OBJ files under data/.
2. Each mesh vertex count is within ±5% of the target.
3. Each mesh has face/vertex Euler ratio in [1.8, 2.2] (subdivision pedigree).

These tests are marked 'slow' because mesh generation takes time.
After the script has been run once and OBJ files committed, re-running
the tests is fast (they only load/check the existing files).
"""

from pathlib import Path

import numpy as np
import pytest

DATA_DIR = Path(__file__).parent.parent / "data"

EXPECTED_MESHES = {
    "Template_1K": 1_000,
    "Template_10K": 10_000,
    "Template_100K": 100_000,
    "DemoFace_1K": 1_000,
    "DemoFace_10K": 10_000,
    "DemoFace_100K": 100_000,
}

TOLERANCE_PCT = 5.0
EULER_RATIO_MIN = 1.8
EULER_RATIO_MAX = 2.2


@pytest.fixture(scope="module")
def trimesh_module():
    """Import trimesh, skip if unavailable."""
    trimesh = pytest.importorskip("trimesh")
    return trimesh


@pytest.mark.parametrize("mesh_name,target_n", list(EXPECTED_MESHES.items()))
def test_profiling_mesh_exists(mesh_name, target_n):
    """Each generated OBJ file must exist on disk."""
    obj_path = DATA_DIR / f"{mesh_name}.obj"
    assert obj_path.exists(), (
        f"Missing profiling mesh: {obj_path}\n"
        f"Run: python scripts/generate_profiling_meshes.py"
    )


@pytest.mark.parametrize("mesh_name,target_n", list(EXPECTED_MESHES.items()))
def test_profiling_mesh_vertex_count(mesh_name, target_n, trimesh_module):
    """Each mesh vertex count must be within ±5% of its target."""
    obj_path = DATA_DIR / f"{mesh_name}.obj"
    if not obj_path.exists():
        pytest.skip(f"OBJ not yet generated: {obj_path}")
    mesh = trimesh_module.load(str(obj_path), process=False)
    n = len(mesh.vertices)
    diff_pct = abs(n - target_n) / target_n * 100
    assert diff_pct <= TOLERANCE_PCT, (
        f"{mesh_name}: {n} verts is {diff_pct:.2f}% from target {target_n} "
        f"(max allowed: {TOLERANCE_PCT}%)"
    )


@pytest.mark.parametrize("mesh_name,target_n", list(EXPECTED_MESHES.items()))
def test_profiling_mesh_euler_ratio(mesh_name, target_n, trimesh_module):
    """Each mesh must have F/V ratio in [1.8, 2.2] (genus-0 triangulated mesh).

    This confirms the subdivision pedigree was used (ADR-004 D1 FIRM):
    a properly subdivided and then decimated mesh preserves the ~2.0 ratio.
    """
    obj_path = DATA_DIR / f"{mesh_name}.obj"
    if not obj_path.exists():
        pytest.skip(f"OBJ not yet generated: {obj_path}")
    mesh = trimesh_module.load(str(obj_path), process=False)
    n_verts = len(mesh.vertices)
    n_faces = len(mesh.faces)
    ratio = n_faces / n_verts
    assert EULER_RATIO_MIN <= ratio <= EULER_RATIO_MAX, (
        f"{mesh_name}: F/V ratio={ratio:.3f} is outside [{EULER_RATIO_MIN}, {EULER_RATIO_MAX}]. "
        f"({n_verts} verts, {n_faces} faces). "
        f"Check that subdivision pedigree was used (not raw-mesh decimation)."
    )
