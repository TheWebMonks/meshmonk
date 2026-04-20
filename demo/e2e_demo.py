"""End-to-end demo: chained rigid → nonrigid → pyramid.

Each stage takes the *previous* stage's output as its input, so the
per-stage plots show the incremental improvement at each step.

Usage (from repo root or demo/):
    cd demo && uv run e2e_demo.py

Inputs:  Template.obj   (floating — the template face)
         demoFace.obj   (target  — the subject face)
Outputs:
    e2e_rigid.obj       — after rigid stage
    e2e_nonrigid.obj    — after nonrigid stage (input = rigid output)
    e2e_pyramid.obj     — after pyramid stage  (input = nonrigid output)
    e2e_rigid.png       — original   → rigid
    e2e_nonrigid.png    — rigid      → nonrigid
    e2e_pyramid.png     — nonrigid   → pyramid
    e2e_compare.png     — side-by-side progression across all stages
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np

import meshmonk
from _common import (
    FLOATING,
    TARGET,
    _view_axes,
    _zoom_bounds,
    load_meshes,
    multi_view_plot,
    quiet_stderr,
    save_obj,
    summarise,
)

HERE = Path(__file__).parent


def _compare_plot(floating_v, target_v, stages: list[dict], out_png: Path) -> None:
    """3x(1+N) grid: rows = planes, cols = original + each stage (cumulative)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    planes = ("front", "side", "top")
    ncols = 1 + len(stages)
    fig, axes = plt.subplots(len(planes), ncols, figsize=(4 * ncols, 4 * len(planes)))

    for row, plane in enumerate(planes):
        ix, iy, xlabel, ylabel = _view_axes(plane)

        ax0 = axes[row, 0]
        ax0.scatter(target_v[:, ix], target_v[:, iy], s=0.3, c="tomato",
                    alpha=0.35, rasterized=True)
        ax0.scatter(floating_v[:, ix], floating_v[:, iy], s=0.4, c="steelblue",
                    alpha=0.7, rasterized=True)
        ax0.set_title(f"original — {plane}")
        ax0.set_xlabel(xlabel); ax0.set_ylabel(ylabel)
        ax0.set_aspect("equal")
        (xlo, xhi), (ylo, yhi) = _zoom_bounds(floating_v, ix, iy)
        ax0.set_xlim(xlo, xhi); ax0.set_ylim(ylo, yhi)

        for col, stage in enumerate(stages, start=1):
            ax = axes[row, col]
            aligned = stage["aligned"]
            ax.scatter(target_v[:, ix], target_v[:, iy], s=0.3, c="tomato",
                       alpha=0.35, rasterized=True)
            ax.scatter(aligned[:, ix], aligned[:, iy], s=0.4, c="seagreen",
                       alpha=0.7, rasterized=True)
            ax.set_title(f"{stage['name']} — {plane}")
            ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
            ax.set_aspect("equal")
            (xlo, xhi), (ylo, yhi) = _zoom_bounds(aligned, ix, iy)
            ax.set_xlim(xlo, xhi); ax.set_ylim(ylo, yhi)

    fig.suptitle("End-to-end progression (cumulative)", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(str(out_png), dpi=120)
    plt.close(fig)
    print(f"  saved compare plot → {out_png}")


def main() -> None:
    meshmonk.set_log_level("warning")

    print("\n=== meshmonk end-to-end demo (chained) ===\n")
    print("Loading meshes...")
    floating, target = load_meshes()
    summarise(floating, "floating (Template)")
    summarise(target, "target  (demoFace)")

    floating_v = np.asarray(floating.vertices, dtype=np.float32)
    target_v_f = np.asarray(target.vertices, dtype=np.float32)
    floating_f = np.asarray(floating.faces, dtype=np.int32)
    target_f = np.asarray(target.faces, dtype=np.int32)

    # Build Pattern-B feature matrices once. We feed each stage's aligned_features
    # forward as the next stage's floating_features.
    feat_float = meshmonk.features_from_vertices(floating_v, floating_f)
    feat_target = meshmonk.features_from_vertices(target_v_f, target_f)
    target_v = target_v_f  # for plotting

    stages: list[dict] = []

    # -- Stage 1: rigid (input = original floating) --
    print("\nStage 1/3 — rigid (SE(3) + scale)...")
    t0 = time.perf_counter()
    with quiet_stderr():
        r_rigid = meshmonk.rigid_register(
            floating_features=feat_float,
            target_features=feat_target,
            floating_faces=floating_f,
            target_faces=target_f,
            num_iterations=40,
            use_scaling=True,
        )
    elapsed = time.perf_counter() - t0
    aligned_rigid = r_rigid.aligned_vertices
    feat_after_rigid = r_rigid.aligned_features
    print(f"  done in {elapsed:.1f}s  iterations={r_rigid.iterations_run}")
    save_obj(aligned_rigid, floating_f, HERE / "e2e_rigid.obj")
    multi_view_plot(
        floating_v, target_v, aligned_rigid, HERE / "e2e_rigid.png",
        title="Stage 1: rigid  (original → rigid)",
    )
    stages.append({"name": "rigid", "aligned": aligned_rigid, "elapsed": elapsed})

    # -- Stage 2: nonrigid (input = rigid output) --
    print("\nStage 2/3 — nonrigid (input = rigid output)...")
    t0 = time.perf_counter()
    with quiet_stderr():
        r_nonrigid = meshmonk.nonrigid_register(
            floating_features=feat_after_rigid,
            target_features=feat_target,
            floating_faces=floating_f,
            target_faces=target_f,
            num_iterations=40,
        )
    elapsed = time.perf_counter() - t0
    aligned_nonrigid = r_nonrigid.aligned_vertices
    feat_after_nonrigid = r_nonrigid.aligned_features
    print(f"  done in {elapsed:.1f}s  mean-inlier={r_nonrigid.final_inlier_weights.mean():.3f}")
    save_obj(aligned_nonrigid, floating_f, HERE / "e2e_nonrigid.obj")
    multi_view_plot(
        aligned_rigid, target_v, aligned_nonrigid, HERE / "e2e_nonrigid.png",
        title="Stage 2: nonrigid  (rigid → nonrigid)",
    )
    stages.append({"name": "nonrigid", "aligned": aligned_nonrigid, "elapsed": elapsed})

    # -- Stage 3: pyramid (input = nonrigid output) --
    print("\nStage 3/3 — pyramid (input = nonrigid output)...")
    t0 = time.perf_counter()
    with quiet_stderr():
        r_pyramid = meshmonk.pyramid_register(
            floating_features=feat_after_nonrigid,
            target_features=feat_target,
            floating_faces=floating_f,
            target_faces=target_f,
            num_iterations=10,
            num_pyramid_layers=3,
        )
    elapsed = time.perf_counter() - t0
    aligned_pyramid = r_pyramid.aligned_vertices
    print(f"  done in {elapsed:.1f}s  per-layer={r_pyramid.per_layer_iterations}  "
          f"mean-inlier={r_pyramid.final_inlier_weights.mean():.3f}")
    save_obj(aligned_pyramid, floating_f, HERE / "e2e_pyramid.obj")
    multi_view_plot(
        aligned_nonrigid, target_v, aligned_pyramid, HERE / "e2e_pyramid.png",
        title="Stage 3: pyramid  (nonrigid → pyramid)",
    )
    stages.append({"name": "pyramid", "aligned": aligned_pyramid, "elapsed": elapsed})

    # -- Cumulative comparison --
    print("\nRendering cumulative comparison plot...")
    _compare_plot(floating_v, target_v, stages, HERE / "e2e_compare.png")

    # Summary: movement relative to original floating
    print("\nSummary (displacement relative to original floating):")
    prev = floating_v
    for stage in stages:
        aligned = stage["aligned"]
        total = np.linalg.norm(aligned - floating_v, axis=1)
        step = np.linalg.norm(aligned - prev, axis=1)
        print(f"  {stage['name']:<9} {stage['elapsed']:5.1f}s   "
              f"total(mean={total.mean():6.2f} max={total.max():6.2f})   "
              f"step(mean={step.mean():6.2f} max={step.max():6.2f})")
        prev = aligned

    print("\nDone.\n")


if __name__ == "__main__":
    main()
