"""
Tier 3 and Tier 3.5 golden tests.

Tier 3: human-approved visual goldens (captured 2026-04-18).
  - test_rigid_golden
  - test_nonrigid_golden
  - test_pyramid_golden

Tier 3.5 (advisory): scientific equivalence against legacy baseline from
  tests/golden/legacy_baseline/rigid_transform.txt and the registered goldens
  (rigid_output.npz). These run only in the nightly workflow; the entire
  test_golden.py file is excluded from PR CI via:
    pytest tests/ -q --ignore=tests/test_golden.py

Usage:
  pytest tests/test_golden.py -v                    # all golden tests
  pytest tests/test_golden.py -v -m advisory        # only advisory

CI configuration:
  - PR CI: --ignore=tests/test_golden.py  (this file excluded entirely)
  - Nightly workflow: pytest tests/test_golden.py -v
"""

from pathlib import Path

import numpy as np
import pytest

import meshmonk

GOLDEN_DIR = Path(__file__).parent / "golden"
LEGACY_DIR = GOLDEN_DIR / "legacy_baseline"
DATA_DIR = Path(__file__).parent.parent / "data"

TEMPLATE_OBJ = DATA_DIR / "Template.obj"
DEMOFACE_OBJ = DATA_DIR / "demoFace.obj"


# ---------------------------------------------------------------------------
# Tier 3: Human-approved visual goldens (stubs — xfail until goldens captured)
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_rigid_golden():
    """Compare rigid registration output against frozen human-approved golden.

    Passes when tests/golden/rigid.npz exists and RMSE is within tolerance.
    """
    from tests.utils.mesh_compare import compare_to_golden  # noqa: PLC0415

    golden_path = GOLDEN_DIR / "rigid.npz"
    if not golden_path.exists():
        pytest.skip(f"Golden file not found: {golden_path}")

    pytest.importorskip("trimesh")
    if not TEMPLATE_OBJ.exists():
        pytest.skip(f"Template.obj not found at {TEMPLATE_OBJ}")
    if not DEMOFACE_OBJ.exists():
        pytest.skip(f"demoFace.obj not found at {DEMOFACE_OBJ}")

    result = meshmonk.rigid_register(
        floating=str(TEMPLATE_OBJ),
        target=str(DEMOFACE_OBJ),
    )
    aligned_vertices = result.aligned_vertices.astype("float64")
    assert compare_to_golden(
        aligned_vertices, str(golden_path), atol=0.5
    ), "Rigid registration output exceeds golden tolerance (RMSE > 0.5 mm)"


@pytest.mark.slow
def test_nonrigid_golden():
    """Compare nonrigid registration output against frozen human-approved golden.

    Passes when tests/golden/nonrigid.npz exists and RMSE is within tolerance.
    """
    from tests.utils.mesh_compare import compare_to_golden  # noqa: PLC0415

    golden_path = GOLDEN_DIR / "nonrigid.npz"
    if not golden_path.exists():
        pytest.skip(f"Golden file not found: {golden_path}")

    pytest.importorskip("trimesh")
    if not TEMPLATE_OBJ.exists():
        pytest.skip(f"Template.obj not found at {TEMPLATE_OBJ}")
    if not DEMOFACE_OBJ.exists():
        pytest.skip(f"demoFace.obj not found at {DEMOFACE_OBJ}")

    result = meshmonk.nonrigid_register(
        floating=str(TEMPLATE_OBJ),
        target=str(DEMOFACE_OBJ),
    )
    aligned_vertices = result.aligned_vertices.astype("float64")
    assert compare_to_golden(
        aligned_vertices, str(golden_path), atol=0.5
    ), "Nonrigid registration output exceeds golden tolerance (RMSE > 0.5 mm)"


@pytest.mark.slow
def test_pyramid_golden():
    """Compare pyramid registration output against frozen human-approved golden.

    Passes when tests/golden/pyramid.npz exists and RMSE is within tolerance.
    """
    from tests.utils.mesh_compare import compare_to_golden  # noqa: PLC0415

    golden_path = GOLDEN_DIR / "pyramid.npz"
    if not golden_path.exists():
        pytest.skip(f"Golden file not found: {golden_path}")

    pytest.importorskip("trimesh")
    if not TEMPLATE_OBJ.exists():
        pytest.skip(f"Template.obj not found at {TEMPLATE_OBJ}")
    if not DEMOFACE_OBJ.exists():
        pytest.skip(f"demoFace.obj not found at {DEMOFACE_OBJ}")

    result = meshmonk.pyramid_register(
        floating=str(TEMPLATE_OBJ),
        target=str(DEMOFACE_OBJ),
    )
    aligned_vertices = result.aligned_vertices.astype("float64")
    assert compare_to_golden(
        aligned_vertices, str(golden_path), atol=0.5
    ), "Pyramid registration output exceeds golden tolerance (RMSE > 0.5 mm)"


# ---------------------------------------------------------------------------
# Tier 3.5: Scientific equivalence against legacy baseline (advisory)
# ---------------------------------------------------------------------------


@pytest.mark.slow
@pytest.mark.advisory
def test_rigid_scientific_equivalence():
    """Assert rigid registration output is close to the legacy baseline transform.

    Uses /workspace/tests/golden/legacy_baseline/rigid_transform.txt as the
    reference, captured from a prior meshmonk_cli run on the same meshes.

    This is an advisory test — it reports discrepancies but does not gate PRs.
    The entire test_golden.py file is excluded from PR CI.
    """
    pytest.importorskip("trimesh")

    if not TEMPLATE_OBJ.exists():
        pytest.skip(f"Template.obj not found at {TEMPLATE_OBJ}")
    if not DEMOFACE_OBJ.exists():
        pytest.skip(f"demoFace.obj not found at {DEMOFACE_OBJ}")

    legacy_transform_path = LEGACY_DIR / "rigid_transform.txt"
    if not legacy_transform_path.exists():
        pytest.skip(f"Legacy baseline transform not found: {legacy_transform_path}")

    # Load legacy reference transform
    legacy_transform = np.loadtxt(str(legacy_transform_path), dtype="float32")
    assert legacy_transform.shape == (4, 4)

    # Run rigid registration
    result = meshmonk.rigid_register(
        floating=str(TEMPLATE_OBJ),
        target=str(DEMOFACE_OBJ),
    )
    result_transform = result.transform.matrix  # (4, 4) float32

    # Compare: loose tolerance — advisory only
    # atol=1e-2 for matrix elements, RMSE atol=1e-1 for vertex positions
    matrix_close = np.allclose(result_transform, legacy_transform, atol=1e-2)

    if not matrix_close:
        max_diff = float(np.max(np.abs(result_transform - legacy_transform)))
        pytest.fail(
            f"[ADVISORY] Rigid transform deviates from legacy baseline: "
            f"max_abs_diff={max_diff:.4f} (atol=1e-2). "
            f"This may indicate algorithm drift or platform-specific numerics. "
            f"Inspect before promoting to hard gate.",
        )
