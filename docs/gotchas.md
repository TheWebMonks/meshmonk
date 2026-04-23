# Common Gotchas

This page covers the most frequent mistakes when using MeshMonk. Each section
explains what the gotcha is, why it matters, and how to avoid it.

---

## Float32 — input arrays are silently coerced

**Gotcha:** Feature arrays, face arrays, and flag arrays are all silently cast to
the required dtype (`float32` or `int32`) by the Python layer. Passing a `float64`
feature matrix does not raise an error — it is silently truncated to `float32`
precision.

**Why it matters:** If your coordinates were computed in float64 and you inspect
the returned `aligned_vertices`, you will see float32 precision. This is usually
fine for mesh registration, but if you concatenate the returned vertices with
float64 data without an explicit cast, NumPy will silently truncate.

**How to avoid:** Cast explicitly before and after registration if downstream code
requires float64:

```python
import numpy as np
aligned = np.asarray(result.aligned_vertices, dtype=np.float64)
```

All result fields (`aligned_features`, `aligned_vertices`, `displacement_field`,
`final_inlier_weights`) are always returned as `float32`.

---

## Face dtype and shape — must be triangles, signed int

**Gotcha:** The C++ core requires face arrays to be `(M, 3)` signed int32 — exactly
3 columns (triangles). Quad meshes and polygon soups are not supported. Passing a
face array with 4 or more columns raises a `DegenerateInput` error.

**Why it matters:** Some mesh loaders produce quad meshes (e.g., from CAD software
or subdivision surfaces). `trimesh.load` triangulates automatically for most formats,
but not all.

**How to avoid:** Ensure your mesh is triangulated before registration. With trimesh:

```python
import trimesh

mesh = trimesh.load("my_mesh.obj")
if not mesh.is_volume:
    mesh = mesh.triangulate()  # splits quads into triangles in-place
```

Face indices from trimesh are `int64` by default. The Python wrapper casts them to
`int32` automatically, so you do not need to convert manually when using Pattern A
(mesh objects or file paths).

---

## Normal direction convention — outward-facing expected

**Gotcha:** The inlier detection step (enabled by default via `inlier_use_orientation=True`)
uses the dot product of corresponding normals to score matches. If one mesh has its
normals flipped relative to the other, nearly all correspondences will be penalised
and registration will fail or produce garbage.

**Why it matters:** The library expects outward-facing normals. OpenMesh computes
normals from face winding order: counter-clockwise winding (viewed from outside)
produces outward normals. If your mesh was exported with clockwise winding, normals
will point inward and orientation-based inlier weights will be near zero.

**How to avoid:** If you see very low inlier weights or a MeshMonk warning about
flipped normals, check the winding order of your mesh. You can either fix the mesh
at export time, or disable the orientation check for a quick diagnostic:

```python
result = meshmonk.rigid_register(
    floating="my_mesh.obj",
    target="target.obj",
    inlier_use_orientation=False,   # disable normal-direction penalty
)
```

When using Pattern A, normals are taken from `mesh.vertex_normals` if they are
non-zero, or recomputed via `compute_normals` (which uses the face winding order).
Pass `compute_normals_flag=True` to force recomputation even when mesh normals exist.

When using Pattern B, if the normal columns (3–5) of `floating_features` or
`target_features` are all-zero, the Python wrapper emits a `UserWarning` and
auto-recomputes normals from the corresponding faces before dispatching to C++.
Pre-compute and supply real normals if you want to silence the warning.

---

## Pattern A vs Pattern B — don't mix them

**Gotcha:** Passing any Pattern B argument (`floating_features`, `target_features`,
`floating_faces`, `target_faces`) alongside a Pattern A argument (`floating`,
`target`) raises a `TypeError`. The two patterns are mutually exclusive.

**Why it matters:** It is easy to accidentally combine them, for example when adding
per-vertex flags to a Pattern A call. Flags can be added to Pattern A calls via the
dedicated `floating_flags` / `target_flags` keyword arguments without switching to
Pattern B.

The full explanation with a mapping table of all call patterns is in the
`rigid_register` docstring (see [Rigid Registration API](api/rigid.md)).

Quick summary:

| Pattern | When to use |
|---|---|
| **A** (`floating=`, `target=`) | You have trimesh objects or file paths. Most common. |
| **B** (`floating_features=`, etc.) | You have pre-built `(N, 6)` feature arrays and want full control. |

---

## All-zero flags — raises DegenerateInput

**Gotcha:** Flags are `float32` values in the range `[0.0, 1.0]`. If *all* flags
for either the floating or target mesh are `0.0`, the C++ core raises
`MeshMonkError` with code `RegistrationError.DegenerateInput` before any
registration work is done.

**Why it matters:** Flags of `0.0` mean "masked out / excluded from registration".
If all vertices are masked, there is nothing to register against.

**How to avoid:** Verify that your flag arrays have at least some non-zero entries
before passing them. By default (Pattern A without explicit flags, or Pattern B
without `floating_flags`/`target_flags`) flags are all-ones, so this only bites
when you supply custom flag arrays.

```python
import numpy as np

flags = np.zeros(N, dtype=np.float32)
flags[active_indices] = 1.0
assert flags.max() > 0, "all flags are zero — registration will raise DegenerateInput"
```

Non-binary values between 0.0 and 1.0 are allowed. The correspondence step
computes a weighted average of target features, then thresholds the resulting
interpolated flag value against `correspondences_flag_threshold` (default ~0.9
for rigid/nonrigid, 0.999 for pyramid). Vertices whose interpolated flag falls
below the threshold are treated as masked at that step.

---

## Pyramid-register MATLAB defaults — flag threshold and viscous/elastic start

**Gotcha:** `pyramid_register` silently applies two non-obvious defaults that differ
from `rigid_register` and `nonrigid_register`:

1. `correspondences_flag_threshold` defaults to `0.999` (not the C++ struct default
   of ~0.9) when not explicitly passed.
2. `transform_num_viscous_iterations_start` and `transform_num_elastic_iterations_start`
   are auto-set to `num_iterations` when not explicitly passed.

**Why it matters:** These match the original MATLAB MeshMonk convention. If you
compare `pyramid_register` output with `nonrigid_register` output using the same
kwargs, the flag threshold and iteration schedule will differ unless you pass the
values explicitly.

**How to avoid:** Pass the values explicitly if you want identical behaviour to
nonrigid:

```python
result = meshmonk.pyramid_register(
    floating="Template.obj",
    target="demoFace.obj",
    correspondences_flag_threshold=0.9,            # override MATLAB default
    transform_num_viscous_iterations_start=10,     # override auto-set
    transform_num_elastic_iterations_start=10,
    num_iterations=100,
)
```
