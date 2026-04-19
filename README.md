# MeshMonk: High-Performance 3D Mesh Registration

> **MATLAB users:** The MATLAB MEX interface has been removed from this repository as part of the Python-first rewrite.
> For MATLAB support, use the university fork at https://github.com/TheWebMonks/meshmonk

MeshMonk is a C++20 library with Python bindings for high-performance 3D mesh registration. It provides rigid, non-rigid, and pyramid-based registration via a clean Python API backed by a battle-tested C++ algorithmic core.

## Overview

MeshMonk v0.1 is a Python-first redesign. The library can be used:

- **As a Python library:** `pip install 'meshmonk[io]'` — then `import meshmonk`
- **As a CLI tool:** `meshmonk rigid Template.obj demoFace.obj --out result.obj`

## Installation

**From PyPI** (when available):

```bash
pip install 'meshmonk[io]'
```

**From source** (current method):

```bash
git clone https://github.com/jsnyde0/meshmonk.git
cd meshmonk
uv pip install --system '.[io]'
```

## Python API Quickstart

```python
import meshmonk, trimesh

floating = trimesh.load("Template.obj")
target   = trimesh.load("demoFace.obj")

result = meshmonk.pyramid_register(floating=floating, target=target)
print(result.aligned_vertices.shape)  # (N, 3)
```

## CLI Usage

```bash
# Rigid registration
meshmonk rigid Template.obj demoFace.obj --out result.obj

# Pyramid (multi-resolution nonrigid) registration
meshmonk pyramid Template.obj demoFace.obj --out result.obj

# View all options for a subcommand
meshmonk rigid --help
meshmonk pyramid --help
```

## CI Matrix

The test suite runs across the following platforms on every pull request:

| OS | Python | Compiler |
|---|---|---|
| Ubuntu 22.04 | 3.10, 3.11, 3.12, 3.13 | gcc-11, clang-15 |
| Ubuntu 24.04 | 3.10, 3.11, 3.12, 3.13 | gcc-13 |
| macOS 13 (x86_64) | 3.10, 3.11, 3.12, 3.13 | AppleClang |
| macOS 14 (arm64) | 3.10, 3.11, 3.12, 3.13 | AppleClang |
| Windows latest | 3.12 | MSVC |

## Repository Structure

- `CMakeLists.txt`: Root CMake file for the project.
- `library/`: Core MeshMonk C++ library (`libmeshmonk_shared`).
- `bindings/`: nanobind Python bindings.
- `meshmonk/`: Python package (`import meshmonk`).
- `vendor/`: Third-party dependencies (OpenMesh, nanoflann, Eigen).
- `tests/`: Python test suite.
- `docs/`: Documentation and architecture decisions.

## Development Setup

```bash
git clone https://github.com/jsnyde0/meshmonk.git
cd meshmonk
uv pip install --system '.[dev,io]'
pytest tests/
```

## Documentation

- [Migrating from MATLAB](docs/migration-from-matlab.md)
- [Architecture decisions](docs/decisions/ADR-001-meshmonk-modernization.md)

## MATLAB Users

The MATLAB MEX interface has been removed in v0.1 as part of the Python-first rewrite (see [ADR-001 D3](docs/decisions/ADR-001-meshmonk-modernization.md)).

For MATLAB support, use the university fork: https://github.com/TheWebMonks/meshmonk
