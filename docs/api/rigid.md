# Rigid Registration

Rigid registration finds the best SE(3) transform — rotation + translation, with
optional uniform scaling — that aligns the floating mesh to the target. Use it when
you expect no surface deformation: pose correction, symmetry alignment, or as a
pre-alignment step before nonrigid or pyramid registration.

It is the fastest of the three registration methods. If the meshes require local
surface deformation to match, use `nonrigid_register` or `pyramid_register` instead.

---

## Worked Example

```python
import meshmonk

# Pattern A: pass file paths directly (requires meshmonk[io])
result = meshmonk.rigid_register(
    floating="Template.obj",
    target="demoFace.obj",
    num_iterations=80,
    use_scaling=False,
)

# Inspect the result
print(result.aligned_vertices.shape)   # (N, 3) float32 — aligned positions
print(result.aligned_features.shape)   # (N, 6) float32 — positions + normals
print(result.transform.matrix)         # (4, 4) float32 — SE(3) homogeneous matrix
print(result.iterations_run)           # int — how many ICP iterations ran

# Save the result
import trimesh, numpy as np
floating = trimesh.load("Template.obj")
registered = trimesh.Trimesh(
    vertices=np.asarray(result.aligned_vertices, dtype=np.float64),
    faces=floating.faces,
    process=False,
)
registered.export("Template_rigid.obj")
```

---

## Return Value — `RigidRegResult`

| Field | Shape / Type | Description |
|---|---|---|
| `aligned_features` | `(N, 6)` float32 | Transformed floating-mesh feature matrix. Columns 0–2 are positions, columns 3–5 are normals. |
| `aligned_vertices` | `(N, 3)` float32 | Shorthand property — `aligned_features[:, :3]`. |
| `transform` | `RigidTransform` | SE(3) rigid transform that maps original floating positions to aligned positions. Access `.matrix` for the `(4, 4)` homogeneous form. |
| `iterations_run` | `int` | Number of ICP iterations actually performed (may be less than `num_iterations` if convergence was early). |

See [Types](types.md) for the full `RigidTransform` API.

---

## Common Pitfalls

See [Gotchas](../gotchas.md) for issues that affect all registration methods, including
float32 requirements, face dtype, and flag semantics.

---

::: meshmonk.rigid_register

::: meshmonk.compute_rigid_transform

::: meshmonk.compute_correspondences

::: meshmonk.compute_inlier_weights
