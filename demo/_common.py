"""Shared helpers for the demo scripts: mesh I/O, quiet stderr, multi-view plots."""

from __future__ import annotations

import contextlib
import os
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
FLOATING = HERE / "Template.obj"
TARGET = HERE / "demoFace.obj"


@contextlib.contextmanager
def quiet_stderr():
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


def summarise(mesh, label: str) -> None:
    v = np.asarray(mesh.vertices)
    lo = v.min(axis=0).round(2)
    hi = v.max(axis=0).round(2)
    print(f"  {label}: {v.shape[0]:,} vertices  bbox {lo} → {hi}")


def load_meshes(floating_path: Path = FLOATING, target_path: Path = TARGET):
    import trimesh
    floating = trimesh.load(str(floating_path))
    target = trimesh.load(str(target_path))
    return floating, target


def save_obj(vertices: np.ndarray, faces: np.ndarray, path: Path) -> None:
    import trimesh
    out = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
    out.export(str(path))


def _view_axes(plane: str):
    """Return (ix, iy, xlabel, ylabel) for a given plane name."""
    if plane == "front":  # XY (looking down -Z)
        return 0, 1, "X", "Y"
    if plane == "side":   # ZY (looking down +X)
        return 2, 1, "Z", "Y"
    if plane == "top":    # XZ (looking down -Y)
        return 0, 2, "X", "Z"
    raise ValueError(f"unknown plane {plane!r}")


def _zoom_bounds(points: np.ndarray, ix: int, iy: int, pad: float = 0.10):
    """Bounds for a 2D projection, padded by `pad` of the max extent."""
    lo = points.min(axis=0)
    hi = points.max(axis=0)
    span = max(hi[ix] - lo[ix], hi[iy] - lo[iy])
    pad_abs = span * pad
    cx = 0.5 * (lo[ix] + hi[ix])
    cy = 0.5 * (lo[iy] + hi[iy])
    half = 0.5 * span + pad_abs
    return (cx - half, cx + half), (cy - half, cy + half)


def multi_view_plot(
    floating_v: np.ndarray,
    target_v: np.ndarray,
    aligned_v: np.ndarray,
    out_png: Path,
    *,
    title: str = "Registration",
    zoom_to: str = "aligned",
):
    """Render a 3x2 grid: front/side/top × (before/after), zoomed on the floating surface.

    `zoom_to` selects which point set drives the viewport bounds. "aligned" zooms
    the after-plots to the aligned floating mesh (the interesting region), and
    the before-plots to the original floating mesh.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    planes = ("front", "side", "top")
    fig, axes = plt.subplots(len(planes), 2, figsize=(10, 12))

    for row, plane in enumerate(planes):
        ix, iy, xlabel, ylabel = _view_axes(plane)

        ax_before = axes[row, 0]
        ax_after = axes[row, 1]

        # Before
        ax_before.scatter(target_v[:, ix], target_v[:, iy], s=0.3, c="tomato",
                          alpha=0.35, label="target", rasterized=True)
        ax_before.scatter(floating_v[:, ix], floating_v[:, iy], s=0.4, c="steelblue",
                          alpha=0.7, label="floating", rasterized=True)
        ax_before.set_title(f"Before — {plane}")
        ax_before.set_xlabel(xlabel); ax_before.set_ylabel(ylabel)
        ax_before.set_aspect("equal")
        (xlo, xhi), (ylo, yhi) = _zoom_bounds(floating_v, ix, iy)
        ax_before.set_xlim(xlo, xhi); ax_before.set_ylim(ylo, yhi)
        if row == 0:
            ax_before.legend(markerscale=12, loc="upper right", fontsize=8)

        # After
        ax_after.scatter(target_v[:, ix], target_v[:, iy], s=0.3, c="tomato",
                         alpha=0.35, label="target", rasterized=True)
        ax_after.scatter(aligned_v[:, ix], aligned_v[:, iy], s=0.4, c="seagreen",
                         alpha=0.7, label="aligned", rasterized=True)
        ax_after.set_title(f"After — {plane}")
        ax_after.set_xlabel(xlabel); ax_after.set_ylabel(ylabel)
        ax_after.set_aspect("equal")
        zoom_pts = aligned_v if zoom_to == "aligned" else floating_v
        (xlo, xhi), (ylo, yhi) = _zoom_bounds(zoom_pts, ix, iy)
        ax_after.set_xlim(xlo, xhi); ax_after.set_ylim(ylo, yhi)
        if row == 0:
            ax_after.legend(markerscale=12, loc="upper right", fontsize=8)

    fig.suptitle(title, fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(str(out_png), dpi=120)
    plt.close(fig)
    print(f"  saved plot → {out_png}")
