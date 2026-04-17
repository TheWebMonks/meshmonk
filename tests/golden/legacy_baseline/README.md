# Legacy Baseline Golden Fixtures

## Schema

All .npz files use the pinned schema from tests/utils/obj_to_npz.py:
- vertices: (N, 3) float64
- faces: (M, 3) int32
- normals: (N, 3) float64 (optional — absent key means not computed)

## Provenance

Source OBJs committed to origin/cli:demo/ from a prior meshmonk_cli run.
Source commit SHA: 0be09c2653cbcf771252070c79342c131cbd3bbc

Equivalent MATLAB demos: test_rigid_registration.m, test_pyramid_registration.m

## NPZ summary (ingested 2026-04-17)

- rigid_output.npz: vertices=(7160, 3), faces=(14050, 3)
- pyramid_output.npz: vertices=(7160, 3), faces=(14050, 3)
- mean per-vertex displacement rigid vs pyramid: 11.2201 mm (> 0.1 mm threshold)

## Tolerance constants

From tests/utils/mesh_compare.py (to be calibrated in v0.1):
- DEFAULT_TRANSFORM_ATOL = 1e-4  (rigid transform components)
- DEFAULT_VERTEX_ATOL_MM = 0.01  (per-vertex position in mm)
