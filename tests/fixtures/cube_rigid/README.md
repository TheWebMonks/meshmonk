# Tier 2 Analytical Fixture: Cube Rigid Transform

## Provenance

Source: `cli/test_data/rigid_transform/` on branch `origin/add-compute-rigid-transform-cli`

The files were cherry-picked from that branch SHA without modification during
bead `meshmonk-modernization-ula.1` (v0.0 reality capture).

## File Manifest

| File | Description |
|---|---|
| `input_mesh.obj` | 8-vertex unit cube OBJ — the floating mesh |
| `input_vertices.txt` | Vertex list in plain text (redundant with OBJ, present upstream) |
| `expected_vertices.txt` | Target vertex positions after the known rigid transform |
| `expected_transform.txt` | 4×4 homogeneous transform matrix (ground truth SE(3)) |
| `corresponding_features.txt` | Per-vertex correspondence features used during `compute_rigid_transform` |
| `inlier_weights.txt` | Per-vertex inlier weights for the rigid fit |

## Purpose

This fixture exercises `compute_rigid_transform` as a **least-squares primitive**,
not end-to-end ICP. Given known correspondences and inlier weights, the function
should recover the expected 4×4 SE(3) transform.

This is a Tier 2 test: grounded in a deterministic geometric fixture rather than
a full registration pipeline. It supplements the Tier 1 synthetic test
(`tests/test_rigid_analytical.py`) which validates end-to-end ICP on a random
point cloud.
