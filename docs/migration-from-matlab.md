# Migrating from MATLAB to Python

MeshMonk v0.1+ replaces the MATLAB MEX interface with a Python API. This guide shows equivalent patterns for common workflows.

## Loading Meshes

**MATLAB** — vertices and faces loaded as separate matrices:

```matlab
% Loaded manually or via custom helpers
template_vertices = ...; % (N, 3) double
template_faces    = ...; % (M, 3) int
```

**Python** — use trimesh (included in `meshmonk[io]`):

```python
import trimesh

mesh = trimesh.load("Template.obj")
vertices = mesh.vertices   # (N, 3) numpy array
faces    = mesh.faces      # (M, 3) numpy array
```

---

## Rigid Registration

**MATLAB**:

```matlab
[aligned_vertices, transform] = meshmonk_rigid( ...
    template_vertices, template_faces, ...
    target_vertices, target_faces, ...
    'numIterations', 80, 'kappa', 12.0);
```

**Python**:

```python
import meshmonk, trimesh

template = trimesh.load("Template.obj")
target   = trimesh.load("demoFace.obj")

result = meshmonk.rigid_register(floating=template, target=target)

aligned_vertices = result.aligned_vertices   # (N, 3) numpy array
transform_matrix = result.transform.matrix  # (4, 4) rigid transform
```

---

## Nonrigid (Viscoelastic) Registration

**MATLAB**:

```matlab
aligned_vertices = meshmonk_nonrigid( ...
    template_vertices, template_faces, ...
    target_vertices, target_faces, ...
    'numIterations', 200);
```

**Python**:

```python
result = meshmonk.nonrigid_register(floating=template, target=target)

aligned_vertices = result.aligned_vertices
displacement     = result.displacement_field  # (N, 3) per-vertex displacement
```

---

## Pyramid (Multi-Resolution Nonrigid) Registration

**MATLAB**:

```matlab
aligned_vertices = meshmonk_pyramid( ...
    template_vertices, template_faces, ...
    target_vertices, target_faces);
```

**Python**:

```python
result = meshmonk.pyramid_register(floating=template, target=target)

aligned_vertices = result.aligned_vertices
displacement     = result.displacement_field  # (N, 3) per-vertex displacement
```

---

## Error Handling (v0.2+ change)

In v0.2, `MeshMonkError` is a dedicated exception class (not just `RuntimeError`).

Catch `meshmonk.MeshMonkError` for intentional library failures. The error message string contains the error type name (e.g. `"DegenerateInput: ..."`).

`RuntimeError` from OpenMesh, Eigen, or other dependencies is no longer caught by `MeshMonkError` — it surfaces directly.

```python
import meshmonk

try:
    result = meshmonk.rigid_register(floating=template, target=target)
except meshmonk.MeshMonkError as e:
    print(f"Registration failed: {e}")
```

---

## Parameter Mapping

| MATLAB name | Python kwarg |
|---|---|
| `numIterations` | `num_iterations` |
| `kappa` | `inlier_kappa` |
| `numNeighbours` | `correspondences_num_neighbours` |
| `usePushPull` | `correspondences_equalize_push_pull` |
| `useScaling` | `use_scaling` |
| `sigma` | `transform_sigma` |
| `numSmoothingNeighbours` | passed inside `NonrigidParams.transform` |
| `numViscousIterations` | `transform_num_viscous_iterations_end` (nonrigid) |
| `numElasticIterations` | `transform_num_elastic_iterations_end` (nonrigid) |

---

## Result Object Fields

| Registration type | Field | Shape | Description |
|---|---|---|---|
| All | `.aligned_vertices` | (N, 3) | Aligned positions |
| All | `.aligned_features` | (N, 6) | Positions + normals |
| Rigid | `.transform.matrix` | (4, 4) | SE(3) rigid transform |
| Nonrigid / Pyramid | `.displacement_field` | (N, 3) | Per-vertex displacement from original |
