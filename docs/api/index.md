# API Overview

MeshMonk exposes two usage patterns:

## Pattern A — convenience functions (recommended)

Pass trimesh `Trimesh` objects directly. MeshMonk extracts vertices and faces automatically.

```python
import trimesh, meshmonk

floating = trimesh.load("floating.obj")
target   = trimesh.load("target.obj")

result = meshmonk.rigid_register(floating=floating, target=target)
```

## Pattern B — raw numpy arrays

Pass `numpy` arrays for vertices and faces explicitly. Useful when working with data
that is not in trimesh format.

```python
import numpy as np, meshmonk

result = meshmonk.rigid_register(
    floating_vertices=np.asarray(..., dtype=np.float32),
    floating_faces=np.asarray(..., dtype=np.int32),
    target_vertices=np.asarray(..., dtype=np.float32),
    target_faces=np.asarray(..., dtype=np.int32),
)
```

---

## Parameter structs

Each registration function accepts either keyword arguments directly or a parameter
struct for batch configuration:

| Struct | Used by |
|--------|---------|
| `RigidParams` | `rigid_register`, `compute_rigid_transform` |
| `NonrigidParams` | `nonrigid_register`, `compute_nonrigid_transform` |
| `PyramidParams` | `pyramid_register` |

---

## Result dataclasses

| Class | Returned by |
|-------|-------------|
| `RigidRegResult` | `rigid_register` |
| `NonrigidRegResult` | `nonrigid_register` |
| `PyramidRegResult` | `pyramid_register` |

All result objects expose `aligned_vertices` (numpy array, shape `(N, 3)`).

---

## Pages

- [Rigid registration](rigid.md)
- [Nonrigid registration](nonrigid.md)
- [Pyramid registration](pyramid.md)
- [Types](types.md)
- [Errors](errors.md)
