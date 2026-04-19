# Quickstart

This guide walks through a complete registration workflow in about 5 minutes.

## Prerequisites

```bash
pip install "meshmonk[io]"
```

`trimesh` (pulled in by `meshmonk[io]`) handles mesh loading and saving.

---

## 1. Load two meshes

```python
import trimesh
import meshmonk

floating = trimesh.load("floating.obj")   # mesh to be registered (source)
target   = trimesh.load("target.obj")     # fixed reference mesh
```

---

## 2. Rigid registration

Rigid registration finds the best SE(3) transform (rotation + translation, optionally
scaling) that aligns `floating` to `target`.

```python
result = meshmonk.rigid_register(
    floating=floating,
    target=target,
    num_iterations=80,
    use_scaling=False,
)

print(result.transform)           # RigidTransform: .rotation (3×3), .translation (3,)
print(result.aligned_vertices)    # numpy array, shape (N, 3)
```

Save the registered mesh:

```python
import numpy as np

registered = trimesh.Trimesh(
    vertices=np.asarray(result.aligned_vertices, dtype=np.float64),
    faces=floating.faces,
    process=False,
)
registered.export("registered_rigid.obj")
```

---

## 3. Pyramid (multi-resolution nonrigid) registration

For dense surface matching, use pyramid registration. It runs nonrigid registration
at multiple mesh resolutions to avoid local minima.

```python
result = meshmonk.pyramid_register(
    floating=floating,
    target=target,
    num_iterations=90,
    num_pyramid_layers=3,
)

print(result.aligned_vertices)    # numpy array, shape (N, 3)
print(result.displacement_field)  # per-vertex displacement, shape (N, 3)
```

Save the result:

```python
registered = trimesh.Trimesh(
    vertices=np.asarray(result.aligned_vertices, dtype=np.float64),
    faces=floating.faces,
    process=False,
)
registered.export("registered_pyramid.obj")
```

---

## 4. Demo data

If you want to try registration with built-in demo meshes, use the CLI:

```bash
meshmonk demo --download   # download Template.obj and demoFace.obj
meshmonk demo pyramid      # run pyramid registration on demo data
```

---

## Next steps

- [API Reference](api/index.md) — all functions, parameters, and return types
- [CLI Reference](cli.md) — full command-line interface documentation
- [Migration from MATLAB](migration-from-matlab.md) — porting MATLAB MeshMonk scripts
