# MeshMonk

MeshMonk is a Python library for 3D mesh registration. It provides rigid (SE(3)), nonrigid
(viscoelastic), and pyramid (multi-resolution nonrigid) registration algorithms, with a
Python API and command-line interface.

> **Note:** These docs cover v0.3 and later. For earlier versions, see the
> [migration guide](migration-from-matlab.md).

---

## Installation

```bash
pip install meshmonk
```

For mesh I/O support (loading and saving OBJ/PLY/STL files):

```bash
pip install "meshmonk[io]"
```

---

## Quick example

```python
import trimesh
import meshmonk

floating = trimesh.load("floating.obj")
target = trimesh.load("target.obj")

result = meshmonk.rigid_register(floating=floating, target=target)

print(result.transform)          # RigidTransform (rotation + translation)
print(result.aligned_vertices)   # numpy array, shape (N, 3)
```

See the [Quickstart](quickstart.md) for a fuller tutorial.

---

## What's next

- [Quickstart](quickstart.md) — load meshes, run rigid and pyramid registration, save results
- [API Reference](api/index.md) — full function signatures and parameter documentation
- [CLI Reference](cli.md) — run registration from the command line
