"""Generate profiling meshes at 1K, 10K, 100K vertices via subdivide-then-downsample.

ADR-004 D1 FIRM: ALL tiers use one loop subdivision pass (1K/10K) or two passes (100K)
before downsampling. Never decimate the raw seed mesh directly.

Run from /workspace:
    python scripts/generate_profiling_meshes.py
Requires: trimesh, fast-simplification, meshmonk (installed with -DMESHMONK_PROFILING not required)
Install extras: pip install "meshmonk[profiling]" or pip install trimesh fast-simplification
"""
from pathlib import Path

import numpy as np
import trimesh

import meshmonk

DATA = Path("/workspace/data")

seeds = [
    ("Template", DATA / "Template.obj"),
    ("DemoFace", DATA / "demoFace.obj"),
]

tiers = [
    ("1K",   1,    1_000),
    ("10K",  1,   10_000),
    ("100K", 2,  100_000),
]


def _remove_non_manifold_faces(mesh: trimesh.Trimesh) -> trimesh.Trimesh:
    """Iteratively remove faces that share a non-manifold edge (shared by >2 faces).

    trimesh's subdivide_loop requires a manifold mesh. Some seed meshes (e.g.
    demoFace.obj) contain a small number of non-manifold edges — this step
    removes them so subdivision can proceed.
    """
    for _ in range(10):  # at most 10 passes to handle cascading non-manifold edges
        edges = mesh.edges_sorted
        unique_edges, inverse, counts = np.unique(
            edges, axis=0, return_inverse=True, return_counts=True
        )
        face_edges = inverse.reshape(-1, 3)
        bad_face_mask = np.any((counts > 2)[face_edges], axis=1)
        if not np.any(bad_face_mask):
            break
        good_faces = mesh.faces[~bad_face_mask]
        mesh = trimesh.Trimesh(vertices=mesh.vertices, faces=good_faces, process=True)
    return mesh


def _simplify_to_target(
    mesh: trimesh.Trimesh, target_n: int, tolerance_pct: float = 5.0
) -> trimesh.Trimesh:
    """Simplify mesh to within tolerance_pct of target_n vertices using trimesh.

    Uses simplify_quadric_decimation with binary search over face_count to
    land within ±tolerance_pct of target_n. The F/V ratio of the result is
    approximately 1.9 (valid for open surfaces).
    """
    # For an open triangulated surface, F ≈ 1.9 * V; use this as initial guess.
    fv_ratio = 1.9
    lo_fc = int(target_n * fv_ratio * 0.8)
    hi_fc = int(target_n * fv_ratio * 1.3)

    best = None
    best_err = float("inf")

    # Binary-search-style iteration over face_count
    for _ in range(20):
        fc = (lo_fc + hi_fc) // 2
        candidate = mesh.simplify_quadric_decimation(face_count=fc)
        n = len(candidate.vertices)
        err = abs(n - target_n) / target_n * 100
        if err < best_err:
            best = candidate
            best_err = err
        if best_err <= tolerance_pct:
            break
        if n < target_n:
            lo_fc = fc + 1
        else:
            hi_fc = fc - 1
        if lo_fc > hi_fc:
            break

    return best


def subdivide_and_decimate(seed_path: Path, n_subdivisions: int, target_n: int):
    """Load seed, loop-subdivide n_subdivisions times, then iteratively decimate to target_n.

    Returns (vertices_xyz, faces) as numpy arrays suitable for OBJ export.
    vertices_xyz: (N, 3) float32
    faces: (M, 3) int32
    """
    mesh = trimesh.load(str(seed_path), process=False)

    # Remove non-manifold faces before subdivision (subdivide_loop requires a
    # manifold mesh — demoFace.obj has a small number of non-manifold edges).
    mesh = _remove_non_manifold_faces(mesh)

    # Loop subdivision returns a NEW Trimesh object
    subdivided = mesh.subdivide_loop(iterations=n_subdivisions)

    # Cast to exact dtypes required by meshmonk
    verts = np.asarray(subdivided.vertices, dtype=np.float32)
    faces = np.asarray(subdivided.faces, dtype=np.int32)  # int32, NOT int64

    # Build (N, 6) feature matrix: positions + normals (required by downsample_mesh)
    normals = np.asarray(meshmonk.compute_normals(verts, faces), dtype=np.float32)
    feat = np.hstack([verts, normals]).astype(np.float32)
    flags = np.ones(len(verts), dtype=np.float32)

    # Iterative decimation until at or below target_n.
    # meshmonk.downsample_mesh(ratio) interprets ratio as fraction of vertices to REMOVE.
    # Correct formula: ratio_to_remove = 1.0 - target_n / current_n, capped at 0.5.
    current_n = feat.shape[0]
    while current_n > target_n:
        ratio = 1.0 - target_n / current_n
        ratio = min(ratio, 0.5)  # never remove more than 50% per step
        feat, faces, flags, _ = meshmonk.downsample_mesh(feat, faces, flags, ratio)
        new_n = feat.shape[0]
        if new_n >= current_n:
            break  # stalled — stop (boundary lock or degenerate mesh)
        current_n = new_n

    # Check whether the meshmonk result meets BOTH the vertex count tolerance AND the
    # Euler F/V ratio requirement [1.8, 2.2]. Open meshes (e.g. Template, demoFace)
    # have many boundary vertices that are locked by OpenMesh's decimater, which can
    # cause the F/V ratio to drop below 1.8 at low vertex counts (e.g. Template_1K
    # ends up at F/V ~1.46 because ~536 of 1000 vertices are on the boundary).
    # In either failure case, fall back to trimesh's quadric simplification on the
    # original subdivided mesh, which does not lock boundary vertices.
    euler_ratio = faces.shape[0] / current_n if current_n > 0 else 0.0
    needs_fallback = current_n > int(target_n * 1.05) or not (1.8 <= euler_ratio <= 2.2)
    if needs_fallback:
        print(
            f"  meshmonk.downsample_mesh result: {current_n} verts, "
            f"F/V={euler_ratio:.2f} — outside target or Euler bounds. "
            f"Falling back to trimesh simplify_quadric_decimation."
        )
        result = _simplify_to_target(subdivided, target_n)
        return (
            np.asarray(result.vertices, dtype=np.float32),
            np.asarray(result.faces, dtype=np.int32),
        )

    return feat[:, :3], faces  # xyz only for OBJ export


for seed_name, seed_path in seeds:
    for tier_name, n_sub, target_n in tiers:
        print(f"Generating {seed_name}_{tier_name} ...")
        verts, faces = subdivide_and_decimate(seed_path, n_sub, target_n)
        out_path = DATA / f"{seed_name}_{tier_name}.obj"
        trimesh.Trimesh(vertices=verts, faces=faces).export(str(out_path))
        print(f"  Wrote {out_path} ({len(verts)} vertices)")

print("Done.")
