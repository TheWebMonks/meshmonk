"""
OBJ -> NPZ converter for MeshMonk golden fixtures.

Pinned output schema:
    vertices: (N, 3) float64  — required, C-contiguous
    faces:    (M, 3) int32    — required, C-contiguous; triangles only
    normals:  (N, 3) float64  — optional, C-contiguous; absent key means 'not computed'

All arrays are guaranteed C-contiguous so that downstream Eigen::Map consumers
(v0.1+, ADR-001 D6 zero-copy intent) can map them without hidden copies.

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

    vertices = np.ascontiguousarray(mesh.points, dtype=np.float64)

    # Collect triangle faces from all cell blocks; reject any non-triangle cells
    triangle_faces = []
    total_cells = 0
    for cell_block in mesh.cells:
        total_cells += len(cell_block.data)
        if cell_block.type == "triangle":
            triangle_faces.append(cell_block.data)
    if not triangle_faces:
        raise ValueError(f"No triangle faces found in {obj_path}")
    faces = np.ascontiguousarray(np.concatenate(triangle_faces, axis=0), dtype=np.int32)
    dropped = total_cells - len(faces)
    if dropped > 0:
        raise ValueError(
            f"non-triangle cells encountered ({dropped} dropped); retriangulate upstream"
        )

    arrays = {"vertices": vertices, "faces": faces}

    # Include normals if present; meshio stores OBJ vertex normals under "obj:vn"
    if mesh.point_data and "obj:vn" in mesh.point_data:
        arrays["normals"] = np.ascontiguousarray(
            mesh.point_data["obj:vn"], dtype=np.float64
        )

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
