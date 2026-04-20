"""End-to-end demo: rigid first, then nonrigid vs pyramid as alternatives.

Pyramid is a multi-resolution variant of the viscoelastic nonrigid
registration — not a refinement applied on top of it. So the realistic
shape of the pipeline is:

    rigid  →  nonrigid       (single-resolution viscoelastic)
       or
    rigid  →  pyramid        (multi-resolution viscoelastic)

This demo runs rigid once, then applies both nonrigid and pyramid to
that rigid output so you can compare the two deformable approaches
side-by-side on equal footing.

Usage (from repo root or demo/):
    cd demo && uv run e2e_demo.py

Inputs:  Template.obj   (floating — the template face)
         demoFace.obj   (target  — the subject face)
Outputs:
    e2e_rigid.obj       — after rigid
    e2e_nonrigid.obj    — rigid output, then nonrigid
    e2e_pyramid.obj     — rigid output, then pyramid
    e2e_rigid.png       — original   → rigid
    e2e_nonrigid.png    — rigid      → nonrigid
    e2e_pyramid.png     — rigid      → pyramid
    e2e_compare.png     — side-by-side: original | rigid | nonrigid | pyramid
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np

import meshmonk
from _common import (
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
    """3x(1+N) grid: rows = planes, cols = original + each stage output."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    planes = ("front", "side", "top")
    ncols = 1 + len(stages)
    fig, axes = plt.subplots(len(planes), ncols, figsize=(4 * ncols, 4 * len(planes)))

    for row, plane in enumerate(planes):
        ix, iy, xlabel, ylabel = _view_axes(plane)

        ax0 = axes[row, 0]
        ax0.scatter(
            target_v[:, ix],
            target_v[:, iy],
            s=0.3,
            c="tomato",
            alpha=0.35,
            rasterized=True,
        )
        ax0.scatter(
            floating_v[:, ix],
            floating_v[:, iy],
            s=0.4,
            c="steelblue",
            alpha=0.7,
            rasterized=True,
        )
        ax0.set_title(f"original — {plane}")
        ax0.set_xlabel(xlabel)
        ax0.set_ylabel(ylabel)
        ax0.set_aspect("equal")
        (xlo, xhi), (ylo, yhi) = _zoom_bounds(floating_v, ix, iy)
        ax0.set_xlim(xlo, xhi)
        ax0.set_ylim(ylo, yhi)

        for col, stage in enumerate(stages, start=1):
            ax = axes[row, col]
            aligned = stage["aligned"]
            ax.scatter(
                target_v[:, ix],
                target_v[:, iy],
                s=0.3,
                c="tomato",
                alpha=0.35,
                rasterized=True,
            )
            ax.scatter(
                aligned[:, ix],
                aligned[:, iy],
                s=0.4,
                c="seagreen",
                alpha=0.7,
                rasterized=True,
            )
            ax.set_title(f"{stage['name']} — {plane}")
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_aspect("equal")
            (xlo, xhi), (ylo, yhi) = _zoom_bounds(aligned, ix, iy)
            ax.set_xlim(xlo, xhi)
            ax.set_ylim(ylo, yhi)

    fig.suptitle("End-to-end: rigid, then nonrigid vs pyramid", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(str(out_png), dpi=120)
    plt.close(fig)
    print(f"  saved compare plot → {out_png}")


def main() -> None:
    meshmonk.set_log_level("warning")

    print("\n=== meshmonk end-to-end demo (rigid → {nonrigid, pyramid}) ===\n")
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
    print("\nStage 1 — rigid (SE(3) + scale)...")
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
        floating_v,
        target_v,
        aligned_rigid,
        HERE / "e2e_rigid.png",
        title="Stage 1: rigid  (original → rigid)",
    )
    stages.append({"name": "rigid", "aligned": aligned_rigid, "elapsed": elapsed})

    # -- Stage 2a: nonrigid on rigid output --
    print("\nStage 2a — nonrigid (input = rigid output)...")
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
    print(
        f"  done in {elapsed:.1f}s  mean-inlier={r_nonrigid.final_inlier_weights.mean():.3f}"
    )
    save_obj(aligned_nonrigid, floating_f, HERE / "e2e_nonrigid.obj")
    multi_view_plot(
        aligned_rigid,
        target_v,
        aligned_nonrigid,
        HERE / "e2e_nonrigid.png",
        title="Stage 2a: nonrigid  (rigid → nonrigid)",
    )
    stages.append({"name": "nonrigid", "aligned": aligned_nonrigid, "elapsed": elapsed})

    # -- Stage 2b: pyramid on rigid output (alternative to nonrigid) --
    print("\nStage 2b — pyramid (input = rigid output, alternative to nonrigid)...")
    t0 = time.perf_counter()
    with quiet_stderr():
        r_pyramid = meshmonk.pyramid_register(
            floating_features=feat_after_rigid,
            target_features=feat_target,
            floating_faces=floating_f,
            target_faces=target_f,
            num_iterations=10,
            num_pyramid_layers=3,
        )
    elapsed = time.perf_counter() - t0
    aligned_pyramid = r_pyramid.aligned_vertices
    print(
        f"  done in {elapsed:.1f}s  per-layer={r_pyramid.per_layer_iterations}  "
        f"mean-inlier={r_pyramid.final_inlier_weights.mean():.3f}"
    )
    save_obj(aligned_pyramid, floating_f, HERE / "e2e_pyramid.obj")
    multi_view_plot(
        aligned_rigid,
        target_v,
        aligned_pyramid,
        HERE / "e2e_pyramid.png",
        title="Stage 2b: pyramid  (rigid → pyramid)",
    )
    stages.append({"name": "pyramid", "aligned": aligned_pyramid, "elapsed": elapsed})

    # -- Cumulative comparison --
    print("\nRendering comparison plot...")
    _compare_plot(floating_v, target_v, stages, HERE / "e2e_compare.png")

    # Summary: rigid is cumulative; nonrigid/pyramid are both measured from rigid.
    print("\nSummary:")
    print(f"  {'stage':<9} {'time':>6}   {'baseline':<8}   mean-move  max-move")
    baselines = {
        "rigid": floating_v,
        "nonrigid": aligned_rigid,
        "pyramid": aligned_rigid,
    }
    labels = {"rigid": "original", "nonrigid": "rigid", "pyramid": "rigid"}
    for stage in stages:
        aligned = stage["aligned"]
        base = baselines[stage["name"]]
        move = np.linalg.norm(aligned - base, axis=1)
        print(
            f"  {stage['name']:<9} {stage['elapsed']:5.1f}s   "
            f"{labels[stage['name']]:<8}   "
            f"{move.mean():6.2f}     {move.max():6.2f}"
        )

    print("\nDone.\n")


if __name__ == "__main__":
    main()
