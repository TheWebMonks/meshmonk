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
| `aligned_vertices` | `ndarray (N, 3)` | Transformed floating mesh vertices. |
| `transform` | `RigidTransform` | Best-fit rigid transform. |
| `iterations_run` | `int` | Number of ICP iterations actually performed. |

### NonrigidRegResult

Returned by `nonrigid_register`.

| Attribute | Type | Description |
|-----------|------|-------------|
| `aligned_vertices` | `ndarray (N, 3)` | Deformed floating mesh vertices. |
| `displacement_field` | `ndarray (N, 3)` | Per-vertex displacement from original position (float32). |
| `final_inlier_weights` | `ndarray (N,)` | Per-vertex inlier weights at convergence. |
| `iterations_run` | `int` | Number of iterations performed. |

### PyramidRegResult

Returned by `pyramid_register`.

| Attribute | Type | Description |
|-----------|------|-------------|
| `aligned_vertices` | `ndarray (N, 3)` | Deformed floating mesh vertices. |
| `displacement_field` | `ndarray (N, 3)` | Per-vertex displacement from original position (float32). |
| `final_inlier_weights` | `ndarray (N,)` | Per-vertex inlier weights at convergence. |
| `per_layer_iterations` | `list[int]` | Iterations run at each pyramid layer. |

---

## Transform type

### RigidTransform

Represents a rigid SE(3) transform (optionally with uniform scaling).

| Attribute | Type | Description |
|-----------|------|-------------|
| `rotation` | `ndarray (3, 3)` | Rotation matrix. |
| `translation` | `ndarray (3,)` | Translation vector. |
| `scale` | `float` | Uniform scale factor (1.0 if `use_scaling=False`). |

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
