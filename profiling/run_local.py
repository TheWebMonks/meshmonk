"""Local profiling runner — wraps run_profile but uses the actual meshmonk repo path.

This shim exists because run_profile.py hardcodes /workspace/data (CI path).
On the developer machine, data lives at <repo>/data/.

Usage:
    uv run python -m profiling.run_local --tiers 100k --modes nonrigid \
        --runs 5 --warmup 1 --seed 42 \
        --out docs/perf/neighbourfinder-breakdown-20260423-<sha>.md
"""

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["EIGEN_DONT_PARALLELIZE"] = "1"

import argparse
import datetime
import subprocess
import sys
from pathlib import Path

import numpy as np

# Resolve the repo root (this file is at <repo>/profiling/run_local.py)
REPO_ROOT = Path(__file__).parent.parent.resolve()
DATA_DIR = REPO_ROOT / "data"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="MeshMonk local profiling driver (dev-machine shim).",
    )
    parser.add_argument(
        "--tiers",
        type=lambda s: [t.strip().lower() for t in s.split(",")],
        default=["100k"],
    )
    parser.add_argument(
        "--modes",
        type=lambda s: [m.strip().lower() for m in s.split(",")],
        default=["nonrigid"],
    )
    parser.add_argument("--runs", type=int, default=5)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out", type=str, default=None, metavar="PATH")
    return parser.parse_args(argv)


def load_mesh_tier_local(tier: str) -> tuple:
    import trimesh
    import meshmonk

    tier_map = {
        "1k": ("Template_1K.obj", "DemoFace_1K.obj"),
        "10k": ("Template_10K.obj", "DemoFace_10K.obj"),
        "100k": ("Template_100K.obj", "DemoFace_100K.obj"),
    }
    tier_key = tier.strip().lower()
    if tier_key not in tier_map:
        raise ValueError(f"Unknown tier {tier!r}")

    float_file, target_file = tier_map[tier_key]
    m_float = trimesh.load(str(DATA_DIR / float_file), process=False)
    m_target = trimesh.load(str(DATA_DIR / target_file), process=False)

    vf = np.asarray(m_float.vertices, dtype=np.float32)
    ff = np.asarray(m_float.faces, dtype=np.int32)
    vt = np.asarray(m_target.vertices, dtype=np.float32)
    ft = np.asarray(m_target.faces, dtype=np.int32)

    feat_f = meshmonk.features_from_vertices(vf, ff)
    feat_t = meshmonk.features_from_vertices(vt, ft)
    flags_f = np.ones(len(vf), dtype=np.float32)
    flags_t = np.ones(len(vt), dtype=np.float32)

    return feat_f, ff, flags_f, feat_t, ft, flags_t


def main():
    import meshmonk

    args = parse_args()

    if not meshmonk.profiling_enabled():
        print(
            "ERROR: meshmonk not built with -DMESHMONK_PROFILING=ON",
            file=sys.stderr,
        )
        sys.exit(2)

    # Import helpers from run_profile (non-path-dependent ones)
    from profiling.run_profile import (
        run_registration,
        collect_stats,
        fit_scaling_exponent,
        generate_report,
    )

    # Calibration
    N_CAL = 1_000_000
    cal_runs = [int(meshmonk.profiling_calibrate(N_CAL)) for _ in range(3)]
    overhead_ns = int(np.median(cal_runs))
    overhead_us = overhead_ns / 1000.0
    print(
        f"ScopedTimer overhead: {cal_runs[0]}, {cal_runs[1]}, {cal_runs[2]} ns/scope"
        f" -> median {overhead_ns} ns/scope ({overhead_us:.3f} us/scope)"
    )

    stats_by_tier_mode = {}

    for tier in args.tiers:
        print(f"\nLoading meshes for tier {tier.upper()} ...")
        feat_f, faces_f, flags_f, feat_t, faces_t, flags_t = load_mesh_tier_local(tier)
        n_verts = feat_f.shape[0]
        print(f"  Floating: {n_verts} vertices, Target: {feat_t.shape[0]} vertices")

        for mode in args.modes:
            print(
                f"  Running {mode} at {tier.upper()} "
                f"({args.warmup} warmup + {args.runs} runs) ..."
            )
            all_data = []
            for rep in range(args.warmup + args.runs):
                meshmonk.profiling_reset()
                run_registration(mode, feat_f, faces_f, flags_f, feat_t, faces_t, flags_t)
                data = meshmonk.profiling_dump()
                all_data.append(data)
                tag = "(warmup)" if rep < args.warmup else ""
                print(f"    rep {rep + 1}/{args.warmup + args.runs} {tag} — {len(data)} labels")
            measured = all_data[args.warmup:]
            stats_by_tier_mode[(tier, mode)] = collect_stats(measured, overhead_us)
            print(f"  {mode} at {tier.upper()}: {len(stats_by_tier_mode[(tier, mode)])} labels collected")

    print("\nFitting scaling exponents ...")
    scaling_exponents = {}
    for tier in args.tiers:
        for mode in args.modes:
            key = (tier, mode)
            if key not in stats_by_tier_mode:
                continue
            for label in stats_by_tier_mode[key]:
                lm_key = (label, mode)
                if lm_key not in scaling_exponents:
                    k = fit_scaling_exponent(stats_by_tier_mode, label, mode)
                    scaling_exponents[lm_key] = k

    print("\nGenerating report ...")
    report_md = generate_report(
        args, overhead_ns, stats_by_tier_mode, scaling_exponents, cal_runs=cal_runs
    )

    if args.out:
        out_path = Path(args.out)
    else:
        try:
            git_short = subprocess.check_output(
                ["git", "rev-parse", "--short", "HEAD"],
                cwd=str(REPO_ROOT),
            ).strip().decode()
        except Exception:
            git_short = "unknown"
        today = datetime.date.today().strftime("%Y%m%d")
        out_dir = REPO_ROOT / "docs" / "perf"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"neighbourfinder-breakdown-{today}-{git_short}.md"

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(report_md, encoding="utf-8")
    print(f"\nReport written to: {out_path}")

    # Sanity: verify the 4 bucket labels appear
    missing = []
    for label in [
        "NeighbourFinder::buffer_alloc",
        "NeighbourFinder::query_setup",
        "NeighbourFinder::tree_query",
        "NeighbourFinder::result_copy",
    ]:
        if label not in report_md:
            missing.append(label)

    if missing:
        print(f"\nWARNING: missing bucket labels in report: {missing}")
        sys.exit(1)
    else:
        print("\nAll 4 bucket labels present in report.")


if __name__ == "__main__":
    main()
