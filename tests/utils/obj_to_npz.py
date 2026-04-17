"""
OBJ -> NPZ converter for MeshMonk golden fixtures.

Pinned output schema:
    vertices: (N, 3) float64  — required
    faces:    (M, 3) int32    — required
    normals:  (N, 3) float64  — optional, absent key means 'not computed'

Usage:
    python -m tests.utils.obj_to_npz <input.obj> <output.npz>
"""

import argparse
from pathlib import Path

import meshio
import numpy as np


def convert(obj_path: "str | Path", npz_path: "str | Path") -> None:
    """Convert an OBJ file to a .npz file using the pinned schema."""
    mesh = meshio.read(str(obj_path))

    vertices = np.asarray(mesh.points, dtype=np.float64)

    # Collect triangle faces from all cell blocks
    triangle_faces = []
    for cell_block in mesh.cells:
        if cell_block.type == "triangle":
            triangle_faces.append(cell_block.data)
    if not triangle_faces:
        raise ValueError(f"No triangle faces found in {obj_path}")
    faces = np.concatenate(triangle_faces, axis=0).astype(np.int32)

    arrays = {"vertices": vertices, "faces": faces}

    # Include normals if present
    if mesh.point_data and "Normals" in mesh.point_data:
        arrays["normals"] = np.asarray(mesh.point_data["Normals"], dtype=np.float64)

    np.savez_compressed(str(npz_path), **arrays)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert OBJ mesh to NPZ using pinned schema."
    )
    parser.add_argument("input_obj", help="Path to input .obj file")
    parser.add_argument("output_npz", help="Path to output .npz file")
    args = parser.parse_args()

    Path(args.output_npz).parent.mkdir(parents=True, exist_ok=True)
    convert(args.input_obj, args.output_npz)
    print(f"Converted {args.input_obj} -> {args.output_npz}")


if __name__ == "__main__":
    main()
