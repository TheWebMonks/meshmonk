"""Demo: nonrigid (viscoelastic) registration.

We run rigid alignment first (because in practice you'd never run nonrigid
on un-aligned meshes), save that intermediate, then run nonrigid on the
rigid output. The per-stage plot shows "rigid → nonrigid" so you can see
what the nonrigid step actually contributes on top of rigid.

Usage:
    cd demo && uv run nonrigid_demo.py
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np

import meshmonk
from _common import (
    load_meshes,
    multi_view_plot,
    quiet_stderr,
    save_obj,
    summarise,
)

HERE = Path(__file__).parent
OUT_RIGID = HERE / "nonrigid_input_rigid.obj"  # rigid pre-alignment (nonrigid input)
OUT_MESH = HERE / "nonrigid_output.obj"
OUT_PNG = HERE / "nonrigid_result.png"


def main() -> None:
    meshmonk.set_log_level("warning")

    print("\n=== meshmonk nonrigid demo (rigid → nonrigid) ===\n")
    print("Loading meshes...")
    floating, target = load_meshes()
    summarise(floating, "floating (Template)")
    summarise(target, "target  (demoFace)")

    floating_v = np.asarray(floating.vertices, dtype=np.float32)
    target_v = np.asarray(target.vertices, dtype=np.float32)
    floating_f = np.asarray(floating.faces, dtype=np.int32)
    target_f = np.asarray(target.faces, dtype=np.int32)

    feat_float = meshmonk.features_from_vertices(floating_v, floating_f)
    feat_target = meshmonk.features_from_vertices(target_v, target_f)

    # --- Step 1: rigid pre-alignment ---
    print("\nStep 1 — rigid pre-alignment (SE(3) + scale)...")
    t0 = time.perf_counter()
    with quiet_stderr():
        rigid = meshmonk.rigid_register(
            floating_features=feat_float,
            target_features=feat_target,
            floating_faces=floating_f,
            target_faces=target_f,
            num_iterations=40,
            use_scaling=True,
        )
    print(
        f"  done in {time.perf_counter() - t0:.1f}s  iterations={rigid.iterations_run}"
    )
    aligned_rigid = rigid.aligned_vertices
    save_obj(aligned_rigid, floating_f, OUT_RIGID)
    print(f"  saved rigid intermediate → {OUT_RIGID}")

    # --- Step 2: nonrigid on rigid output ---
    print("\nStep 2 — nonrigid (viscoelastic) on rigid output...")
    t0 = time.perf_counter()
    with quiet_stderr():
        result = meshmonk.nonrigid_register(
            floating_features=rigid.aligned_features,
            target_features=feat_target,
            floating_faces=floating_f,
            target_faces=target_f,
            num_iterations=40,
        )
    elapsed = time.perf_counter() - t0
    print(f"  done in {elapsed:.1f}s  iterations={result.iterations_run}")
    print(f"  mean inlier weight: {result.final_inlier_weights.mean():.3f}")

    aligned_v = result.aligned_vertices
    displ = np.linalg.norm(result.displacement_field, axis=1)
    print(f"  displacement:       mean={displ.mean():.3f}  max={displ.max():.3f}")

    step_move = np.linalg.norm(aligned_v - aligned_rigid, axis=1)
    print(
        f"  step-move (rigid → nonrigid): mean={step_move.mean():.3f}  max={step_move.max():.3f}"
    )

    print("\nSaving output mesh...")
    save_obj(aligned_v, floating_f, OUT_MESH)
    print(f"  saved → {OUT_MESH}")

    print("\nRendering multi-view plot (rigid → nonrigid)...")
    multi_view_plot(
        aligned_rigid,
        target_v,
        aligned_v,
        OUT_PNG,
        title="Nonrigid registration (input: rigid-aligned)",
    )

    print("\nDone.\n")


if __name__ == "__main__":
    main()
