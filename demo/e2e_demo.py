"""End-to-end demo: runs rigid → nonrigid → pyramid and compares results.

Usage (from repo root or demo/):
    cd demo && uv run e2e_demo.py

Inputs:  Template.obj   (floating — the template face)
         demoFace.obj   (target  — the subject face)
Outputs:
    e2e_rigid.obj       — rigid-only aligned template
    e2e_nonrigid.obj    — nonrigid aligned template
    e2e_pyramid.obj     — pyramid aligned template
    e2e_rigid.png       — front/side/top multi-view (rigid)
    e2e_nonrigid.png    — front/side/top multi-view (nonrigid)
    e2e_pyramid.png     — front/side/top multi-view (pyramid)
    e2e_compare.png     — side-by-side comparison across methods
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np

import meshmonk
from _common import (
    FLOATING,
    TARGET,
    load_meshes,
    multi_view_plot,
    quiet_stderr,
    save_obj,
    summarise,
    _view_axes,
    _zoom_bounds,
)

HERE = Path(__file__).parent


def _compare_plot(floating_v, target_v, results: dict, out_png: Path) -> None:
    """3x(1+N) grid: rows = planes (front/side/top), cols = before + each method."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    planes = ("front", "side", "top")
    methods = list(results.keys())
    ncols = 1 + len(methods)
    fig, axes = plt.subplots(len(planes), ncols, figsize=(4 * ncols, 4 * len(planes)))

    for row, plane in enumerate(planes):
        ix, iy, xlabel, ylabel = _view_axes(plane)

        ax0 = axes[row, 0]
        ax0.scatter(target_v[:, ix], target_v[:, iy], s=0.3, c="tomato",
                    alpha=0.35, rasterized=True)
        ax0.scatter(floating_v[:, ix], floating_v[:, iy], s=0.4, c="steelblue",
                    alpha=0.7, rasterized=True)
        ax0.set_title(f"Before — {plane}")
        ax0.set_xlabel(xlabel); ax0.set_ylabel(ylabel)
        ax0.set_aspect("equal")
        (xlo, xhi), (ylo, yhi) = _zoom_bounds(floating_v, ix, iy)
        ax0.set_xlim(xlo, xhi); ax0.set_ylim(ylo, yhi)

        for col, name in enumerate(methods, start=1):
            ax = axes[row, col]
            aligned = results[name]["aligned"]
            ax.scatter(target_v[:, ix], target_v[:, iy], s=0.3, c="tomato",
                       alpha=0.35, rasterized=True)
            ax.scatter(aligned[:, ix], aligned[:, iy], s=0.4, c="seagreen",
                       alpha=0.7, rasterized=True)
            ax.set_title(f"{name} — {plane}")
            ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
            ax.set_aspect("equal")
            (xlo, xhi), (ylo, yhi) = _zoom_bounds(aligned, ix, iy)
            ax.set_xlim(xlo, xhi); ax.set_ylim(ylo, yhi)

    fig.suptitle("End-to-end comparison", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(str(out_png), dpi=120)
    plt.close(fig)
    print(f"  saved compare plot → {out_png}")


def _run_stage(name, fn, floating_v, target_v, floating_f, out_mesh, out_png, title):
    t0 = time.perf_counter()
    with quiet_stderr():
        result = fn()
    elapsed = time.perf_counter() - t0
    aligned_v = result.aligned_vertices
    print(f"  {name}: {elapsed:.1f}s  → {aligned_v.shape[0]:,} verts")
    save_obj(aligned_v, floating_f, out_mesh)
    multi_view_plot(floating_v, target_v, aligned_v, out_png, title=title)
    return {"aligned": aligned_v, "elapsed": elapsed, "result": result}


def main() -> None:
    meshmonk.set_log_level("warning")

    print("\n=== meshmonk end-to-end demo ===\n")
    print("Loading meshes...")
    floating, target = load_meshes()
    summarise(floating, "floating (Template)")
    summarise(target, "target  (demoFace)")

    floating_v = np.asarray(floating.vertices)
    target_v = np.asarray(target.vertices)
    floating_f = np.asarray(floating.faces)

    print("\nRunning all three stages...")
    results: dict[str, dict] = {}

    results["rigid"] = _run_stage(
        "rigid",
        lambda: meshmonk.rigid_register(
            floating=FLOATING, target=TARGET,
            num_iterations=40, use_scaling=True,
        ),
        floating_v, target_v, floating_f,
        HERE / "e2e_rigid.obj", HERE / "e2e_rigid.png",
        "Rigid (SE(3) + scale)",
    )

    results["nonrigid"] = _run_stage(
        "nonrigid",
        lambda: meshmonk.nonrigid_register(
            floating=FLOATING, target=TARGET,
            rigid_params={"num_iterations": 40, "use_scaling": True},
            num_iterations=40,
        ),
        floating_v, target_v, floating_f,
        HERE / "e2e_nonrigid.obj", HERE / "e2e_nonrigid.png",
        "Nonrigid (viscoelastic)",
    )

    results["pyramid"] = _run_stage(
        "pyramid",
        lambda: meshmonk.pyramid_register(
            floating=FLOATING, target=TARGET,
            rigid_params={"num_iterations": 40, "use_scaling": True},
            num_iterations=10, num_pyramid_layers=3,
        ),
        floating_v, target_v, floating_f,
        HERE / "e2e_pyramid.obj", HERE / "e2e_pyramid.png",
        "Pyramid (3 layers)",
    )

    print("\nRendering comparison plot...")
    _compare_plot(floating_v, target_v, results, HERE / "e2e_compare.png")

    # Summary
    print("\nSummary:")
    for name, r in results.items():
        aligned = r["aligned"]
        diff = np.linalg.norm(aligned - floating_v, axis=1)
        print(f"  {name:<8} {r['elapsed']:5.1f}s   mean-move={diff.mean():6.2f}  max-move={diff.max():6.2f}")

    print("\nDone.\n")


if __name__ == "__main__":
    main()
