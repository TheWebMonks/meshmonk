"""End-to-end demo: pyramid registration of demo faces via the Python API.

Usage (from repo root or demo/):
    cd demo && uv run e2e_demo.py

Inputs:  Template.obj   (floating — the template face)
         demoFace.obj   (target  — the subject face)
Output:  e2e_output.obj (aligned template written to disk)
         e2e_result.png (before/after scatter plot)
"""

import contextlib
import os
import time
from pathlib import Path

import numpy as np

import meshmonk

HERE = Path(__file__).parent
FLOATING = HERE / "Template.obj"
TARGET = HERE / "demoFace.obj"
OUT_MESH = HERE / "e2e_output.obj"
OUT_PNG = HERE / "e2e_result.png"


@contextlib.contextmanager
def _quiet_stderr():
    """Redirect C-level stderr to /dev/null (silences OpenMesh topology warnings)."""
    devnull = os.open(os.devnull, os.O_WRONLY)
    old = os.dup(2)
    os.dup2(devnull, 2)
    try:
        yield
    finally:
        os.dup2(old, 2)
        os.close(devnull)
        os.close(old)


def _summarise(mesh, label):
    v = np.asarray(mesh.vertices)
    print(f"  {label}: {v.shape[0]:,} vertices  bbox {v.min(axis=0).round(2)} → {v.max(axis=0).round(2)}")


def _save_obj(vertices, faces, path):
    import trimesh
    out = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
    out.export(str(path))


def _plot(floating_v, target_v, aligned_v, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].set_title("Before registration")
    axes[0].scatter(floating_v[:, 0], floating_v[:, 1], s=0.3, c="steelblue", alpha=0.5, label="floating", rasterized=True)
    axes[0].scatter(target_v[:, 0], target_v[:, 1], s=0.3, c="tomato", alpha=0.5, label="target", rasterized=True)
    axes[0].legend(markerscale=10)
    axes[0].set_aspect("equal")

    axes[1].set_title("After pyramid registration")
    axes[1].scatter(aligned_v[:, 0], aligned_v[:, 1], s=0.3, c="steelblue", alpha=0.5, label="aligned", rasterized=True)
    axes[1].scatter(target_v[:, 0], target_v[:, 1], s=0.3, c="tomato", alpha=0.5, label="target", rasterized=True)
    axes[1].legend(markerscale=10)
    axes[1].set_aspect("equal")

    fig.tight_layout()
    fig.savefig(str(out_png), dpi=120)
    print(f"  saved plot → {out_png}")


def main():
    meshmonk.set_log_level("warning")

    print("\n=== meshmonk e2e demo ===\n")
    print("Loading meshes...")
    import trimesh
    floating = trimesh.load(str(FLOATING))
    target = trimesh.load(str(TARGET))
    _summarise(floating, "floating (Template)")
    _summarise(target,   "target  (demoFace)")

    floating_v = np.asarray(floating.vertices)
    target_v   = np.asarray(target.vertices)
    floating_f = np.asarray(floating.faces)

    print("\nRunning pyramid_register (with rigid pre-alignment)...")
    t0 = time.perf_counter()

    with _quiet_stderr():
        result = meshmonk.pyramid_register(
            floating=FLOATING,
            target=TARGET,
            rigid_params={},
            num_iterations=10,
            num_pyramid_layers=3,
        )

    elapsed = time.perf_counter() - t0
    print(f"  done in {elapsed:.1f}s")
    print(f"  per-layer iterations: {result.per_layer_iterations}")
    print(f"  mean inlier weight:   {result.final_inlier_weights.mean():.3f}")

    aligned_v  = result.aligned_vertices
    displ_mag  = np.linalg.norm(result.displacement_field, axis=1)
    print(f"  displacement field:   mean={displ_mag.mean():.3f}  max={displ_mag.max():.3f}")

    print("\nSaving output mesh...")
    _save_obj(aligned_v, floating_f, OUT_MESH)
    print(f"  saved → {OUT_MESH}")

    print("\nRendering before/after plot...")
    _plot(floating_v, target_v, aligned_v, OUT_PNG)

    print("\nDone.\n")


if __name__ == "__main__":
    main()
