# Types

MeshMonk defines several C++ types exposed to Python via nanobind. Because they are
C-extension types, mkdocstrings may not auto-generate their documentation. Manual
descriptions are provided here.

---

## Result types

### RigidRegResult

Returned by `rigid_register`.

| Attribute | Type | Description |
|-----------|------|-------------|
| `aligned_features` | `ndarray (N, 6)` | Transformed floating mesh features (positions + normals). Primary result field. |
| `aligned_vertices` | `ndarray (N, 3)` | Transformed floating mesh vertices (derived from `aligned_features`). |
| `transform` | `RigidTransform` | Best-fit rigid transform. |
| `iterations_run` | `int` | Number of ICP iterations actually performed. |

### NonrigidRegResult

Returned by `nonrigid_register`.

| Attribute | Type | Description |
|-----------|------|-------------|
| `aligned_features` | `ndarray (N, 6)` | Deformed floating mesh features (positions + normals). Primary result field. |
| `aligned_vertices` | `ndarray (N, 3)` | Deformed floating mesh vertices (derived from `aligned_features`). |
| `displacement_field` | `ndarray (N, 3)` | Per-vertex displacement from original position (float32). |
| `final_inlier_weights` | `ndarray (N,)` | Per-vertex inlier weights at convergence. |
| `iterations_run` | `int` | Number of iterations performed. |

### PyramidRegResult

Returned by `pyramid_register`.

| Attribute | Type | Description |
|-----------|------|-------------|
| `aligned_features` | `ndarray (N, 6)` | Deformed floating mesh features (positions + normals). Primary result field. |
| `aligned_vertices` | `ndarray (N, 3)` | Deformed floating mesh vertices (derived from `aligned_features`). |
| `displacement_field` | `ndarray (N, 3)` | Per-vertex displacement from original position (float32). |
| `final_inlier_weights` | `ndarray (N,)` | Per-vertex inlier weights at convergence. |
| `per_layer_iterations` | `list[int]` | Iterations run at each pyramid layer. |

---

## Transform type

### RigidTransform

Represents a rigid SE(3) transform (optionally with uniform scaling).

| Attribute / Method | Type | Description |
|-----------|------|-------------|
| `.matrix` | `ndarray (4, 4)` float32 | Homogeneous transformation matrix. |
| `.compose(other)` | `RigidTransform` | Compose this transform with another (returns new transform). |
| `.apply(features)` | `ndarray (N, 6)` | Apply transform to a feature array (positions + normals). |
| `.inverse()` | `RigidTransform` | Compute the inverse transform. |

---

## Parameter structs

### RigidParams

Configuration for `rigid_register` / `compute_rigid_transform`.

Key fields: `num_iterations`, `use_scaling`, `inlier_kappa`,
`correspondences_num_neighbours`.

### NonrigidParams

Configuration for `nonrigid_register` / `compute_nonrigid_transform`.

Key fields: `num_iterations`, `inlier_kappa`, `correspondences_num_neighbours`,
`viscous_iterations`, `elastic_iterations`.

### PyramidParams

Configuration for `pyramid_register`.

Key fields: `num_iterations`, `num_pyramid_layers`, `inlier_kappa`,
`correspondences_num_neighbours`.
