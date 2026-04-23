# Nonrigid Registration

Nonrigid (viscoelastic) registration deforms the floating mesh vertex-by-vertex to
match the target surface. It models the floating mesh as a viscoelastic material:
viscous smoothing reduces high-frequency noise early in the solve, elastic regularisation
preserves shape integrity later. Use it when the two meshes share topology but differ
in surface shape — for example, different facial expressions or body poses.

For large shape differences or meshes that are not pre-aligned, prefer `pyramid_register`,
which runs nonrigid registration at multiple resolutions and is more robust to local minima.
Pass `rigid_params={}` to run rigid pre-alignment automatically before the nonrigid solve.

---

## Worked Example

```python
import meshmonk

# With rigid pre-alignment (recommended when meshes are not already roughly aligned)
result = meshmonk.nonrigid_register(
    floating="Template.obj",
    target="demoFace.obj",
    rigid_params={"num_iterations": 40},   # rigid pre-align, then go nonrigid
    num_iterations=100,
    transform_sigma=3.0,
)

# Inspect the result
print(result.aligned_vertices.shape)       # (N, 3) float32 — deformed positions
print(result.displacement_field.shape)     # (N, 3) float32 — per-vertex displacement
print(result.final_inlier_weights.shape)   # (N,)   float32 — convergence quality
print(result.iterations_run)              # int

# Vertices with low inlier weight had poor correspondences at convergence
import numpy as np
low_quality = np.where(result.final_inlier_weights < 0.5)[0]
print(f"{len(low_quality)} vertices had low inlier weight")
```

---

## Return Value — `NonrigidRegResult`

| Field | Shape / Type | Description |
|---|---|---|
| `aligned_features` | `(N, 6)` float32 | Deformed floating-mesh feature matrix. Columns 0–2 are positions, columns 3–5 are normals. |
| `aligned_vertices` | `(N, 3)` float32 | Shorthand property — `aligned_features[:, :3]`. |
| `displacement_field` | `(N, 3)` float32 | Per-vertex displacement from the original (pre-registration) floating positions. Add to original vertices to recover aligned positions. |
| `final_inlier_weights` | `(N,)` float32 | Per-vertex inlier weights at convergence. Values near 1.0 indicate reliable correspondences; values near 0.0 indicate outlier vertices. |
| `iterations_run` | `int` | Number of nonrigid iterations performed. |

---

## Common Pitfalls

See [Gotchas](../gotchas.md) for issues that affect all registration methods, including
float32 requirements, face dtype, flag semantics, and what `final_inlier_weights` near
zero means in practice.

---

::: meshmonk.nonrigid_register

::: meshmonk.compute_nonrigid_transform
