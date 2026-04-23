# Pyramid Registration

Pyramid (multi-resolution nonrigid) registration runs the nonrigid viscoelastic solver
at successively finer mesh resolutions, from a coarse downsampled level down to the full
mesh. This coarse-to-fine strategy avoids local minima and handles larger initial
misalignments than flat nonrigid registration. It is the recommended method for
production use with facial or anatomical surface data.

Use `nonrigid_register` instead when you need finer control over the smoothing schedule
or when the meshes are already closely aligned and a single-resolution solve is sufficient.

---

## Worked Example

```python
import meshmonk

result = meshmonk.pyramid_register(
    floating="Template.obj",
    target="demoFace.obj",
    rigid_params={},           # run rigid pre-alignment with defaults
    num_iterations=90,
    num_pyramid_layers=3,
)

# Inspect the result
print(result.aligned_vertices.shape)       # (N, 3) float32 — final deformed positions
print(result.displacement_field.shape)     # (N, 3) float32 — total per-vertex displacement
print(result.final_inlier_weights.shape)   # (N,)   float32 — quality at convergence
print(result.per_layer_iterations)         # list[int] — iterations run at each layer

# Save the registered mesh
import trimesh, numpy as np
floating = trimesh.load("Template.obj")
registered = trimesh.Trimesh(
    vertices=np.asarray(result.aligned_vertices, dtype=np.float64),
    faces=floating.faces,
    process=False,
)
registered.export("Template_pyramid.obj")
```

---

## Return Value — `PyramidRegResult`

| Field | Shape / Type | Description |
|---|---|---|
| `aligned_features` | `(N, 6)` float32 | Deformed floating-mesh feature matrix at full resolution. Columns 0–2 are positions, columns 3–5 are normals. |
| `aligned_vertices` | `(N, 3)` float32 | Shorthand property — `aligned_features[:, :3]`. |
| `displacement_field` | `(N, 3)` float32 | Total per-vertex displacement accumulated across all pyramid layers. Add to original floating vertices to recover aligned positions. |
| `final_inlier_weights` | `(N,)` float32 | Per-vertex inlier weights at the final (finest) pyramid layer. Values near 1.0 indicate reliable correspondences. |
| `per_layer_iterations` | `list[int]` | Number of nonrigid iterations performed at each pyramid layer, from coarsest to finest. Length equals `num_pyramid_layers`. |

---

## Common Pitfalls

See [Gotchas](../gotchas.md) for issues that affect all registration methods, including
float32 requirements, face dtype, and the MATLAB-convention defaults that `pyramid_register`
applies automatically (flag threshold and viscous/elastic start iterations).

---

::: meshmonk.pyramid_register
