"""Demo: rigid-only registration (SE(3) with optional scaling).

Usage:
    cd demo && uv run rigid_demo.py
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
)

HERE = Path(__file__).parent
OUT_MESH = HERE / "rigid_output.obj"
OUT_PNG = HERE / "rigid_result.png"


def main() -> None:
    meshmonk.set_log_level("warning")

    print("\n=== meshmonk rigid demo ===\n")
    print("Loading meshes...")
    floating, target = load_meshes()
    summarise(floating, "floating (Template)")
    summarise(target, "target  (demoFace)")

    floating_v = np.asarray(floating.vertices)
    target_v = np.asarray(target.vertices)
    floating_f = np.asarray(floating.faces)

    print("\nRunning rigid_register (SE(3) + scaling)...")
    t0 = time.perf_counter()
    with quiet_stderr():
        result = meshmonk.rigid_register(
            floating=FLOATING,
            target=TARGET,
            num_iterations=40,
            use_scaling=True,
        )
    elapsed = time.perf_counter() - t0
    print(f"  done in {elapsed:.1f}s")
    print(f"  iterations run:    {result.iterations_run}")

    aligned_v = result.aligned_vertices
    T = np.asarray(result.transform.matrix)
    print(f"  transform matrix:\n{T.round(3)}")

    print("\nSaving output mesh...")
    save_obj(aligned_v, floating_f, OUT_MESH)
    print(f"  saved → {OUT_MESH}")

    print("\nRendering multi-view plot...")
    multi_view_plot(
        floating_v,
        target_v,
        aligned_v,
        OUT_PNG,
        title="Rigid registration (SE(3) + scale)",
    )

    print("\nDone.\n")


if __name__ == "__main__":
    main()
