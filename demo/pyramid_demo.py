"""Demo: pyramid (multi-resolution) nonrigid registration.

Usage:
    cd demo && uv run pyramid_demo.py
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
OUT_MESH = HERE / "pyramid_output.obj"
OUT_PNG = HERE / "pyramid_result.png"


def main() -> None:
    meshmonk.set_log_level("warning")

    print("\n=== meshmonk pyramid demo ===\n")
    print("Loading meshes...")
    floating, target = load_meshes()
    summarise(floating, "floating (Template)")
    summarise(target, "target  (demoFace)")

    floating_v = np.asarray(floating.vertices)
    target_v = np.asarray(target.vertices)
    floating_f = np.asarray(floating.faces)

    print("\nRunning pyramid_register (3-layer, rigid pre-alignment)...")
    t0 = time.perf_counter()
    with quiet_stderr():
        result = meshmonk.pyramid_register(
            floating=FLOATING,
            target=TARGET,
            rigid_params={"num_iterations": 40, "use_scaling": True},
            num_iterations=10,
            num_pyramid_layers=3,
        )
    elapsed = time.perf_counter() - t0
    print(f"  done in {elapsed:.1f}s")
    print(f"  per-layer iterations: {result.per_layer_iterations}")
    print(f"  mean inlier weight:   {result.final_inlier_weights.mean():.3f}")

    aligned_v = result.aligned_vertices
    displ = np.linalg.norm(result.displacement_field, axis=1)
    print(f"  displacement:         mean={displ.mean():.3f}  max={displ.max():.3f}")

    print("\nSaving output mesh...")
    save_obj(aligned_v, floating_f, OUT_MESH)
    print(f"  saved → {OUT_MESH}")

    print("\nRendering multi-view plot...")
    multi_view_plot(
        floating_v, target_v, aligned_v, OUT_PNG,
        title="Pyramid registration (3 layers)",
    )

    print("\nDone.\n")


if __name__ == "__main__":
    main()
