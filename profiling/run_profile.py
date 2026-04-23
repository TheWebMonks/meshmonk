"""MeshMonk profiling driver.

Run per-label hotspot profiling across 1K/10K/100K mesh tiers and
rigid/nonrigid/pyramid registration modes. Produces a markdown report.

Usage::

    python -m profiling.run_profile \\
        --tiers 1k,10k,100k \\
        --modes rigid,nonrigid,pyramid \\
        --runs 5 \\
        --warmup 1 \\
        --seed 42

See docs/decisions/ADR-004-profiling.md and
history/2026-04-23-profiling-design.md for design context.
"""

# Top-level imports (module level — no meshmonk import yet)
import argparse
import datetime
import hashlib
import os
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np


# meshmonk is imported INSIDE main() after env vars are set — see Step 0.


def parse_args(argv=None) -> argparse.Namespace:
    """Parse command-line arguments.

    Parameters
    ----------
    argv:
        Argument list (defaults to sys.argv[1:]).

    Returns
    -------
    argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description="MeshMonk profiling driver — produces hotspot report.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--tiers",
        type=lambda s: [t.strip().lower() for t in s.split(",")],
        default=["1k", "10k", "100k"],
        metavar="TIERS",
        help="Comma-separated mesh tiers: 1k,10k,100k (default: all three)",
    )
    parser.add_argument(
        "--modes",
        type=lambda s: [m.strip().lower() for m in s.split(",")],
        default=["rigid", "nonrigid", "pyramid"],
        metavar="MODES",
        help="Comma-separated modes: rigid,nonrigid,pyramid (default: all three)",
    )
    parser.add_argument(
        "--runs",
        type=int,
        default=5,
        help="Number of measured runs per (tier, mode) (default: 5)",
    )
    parser.add_argument(
        "--warmup",
        type=int,
        default=1,
        help="Number of warmup runs discarded per (tier, mode) (default: 1)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed (default: 42; reserved for future use)",
    )
    parser.add_argument(
        "--out",
        type=str,
        default=None,
        metavar="PATH",
        help=(
            "Output file path. If not given, auto-derived as "
            "docs/perf/hotspot-profile-YYYYMMDD-<shortsha>.md"
        ),
    )
    return parser.parse_args(argv)


def sha256_file(path: Path) -> str:
    """Return the SHA-256 hex digest of a file.

    Parameters
    ----------
    path:
        Path to the file.

    Returns
    -------
    str
        64-character hex string.
    """
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def load_mesh_tier(tier: str) -> tuple:
    """Load the floating and target meshes for a given tier.

    Parameters
    ----------
    tier:
        One of '1k', '10k', '100k' (case-insensitive).

    Returns
    -------
    tuple
        (feat_float, faces_float, flags_float, feat_target, faces_target, flags_target)
        feat_*:  (N, 6) float32  — positions + normals
        faces_*: (M, 3) int32
        flags_*: (N,)   float32  — all ones
    """
    import trimesh
    import meshmonk

    DATA = Path("/workspace/data")
    tier_map = {
        "1k": ("Template_1K.obj", "DemoFace_1K.obj"),
        "10k": ("Template_10K.obj", "DemoFace_10K.obj"),
        "100k": ("Template_100K.obj", "DemoFace_100K.obj"),
    }
    tier_key = tier.strip().lower()
    if tier_key not in tier_map:
        raise ValueError(f"Unknown tier {tier!r}. Valid: {list(tier_map.keys())}")

    float_file, target_file = tier_map[tier_key]
    m_float = trimesh.load(str(DATA / float_file), process=False)
    m_target = trimesh.load(str(DATA / target_file), process=False)

    vf = np.asarray(m_float.vertices, dtype=np.float32)
    ff = np.asarray(m_float.faces, dtype=np.int32)
    vt = np.asarray(m_target.vertices, dtype=np.float32)
    ft = np.asarray(m_target.faces, dtype=np.int32)

    feat_f = meshmonk.features_from_vertices(vf, ff)
    feat_t = meshmonk.features_from_vertices(vt, ft)
    flags_f = np.ones(len(vf), dtype=np.float32)
    flags_t = np.ones(len(vt), dtype=np.float32)

    return feat_f, ff, flags_f, feat_t, ft, flags_t


def run_registration(
    mode: str,
    feat_f,
    faces_f,
    flags_f,
    feat_t,
    faces_t,
    flags_t,
) -> None:
    """Run a registration and discard the result.

    Profiling data accumulates in the C++ global accumulator. Call
    meshmonk.profiling_reset() before and meshmonk.profiling_dump() after.

    Parameters
    ----------
    mode:
        One of 'rigid', 'nonrigid', 'pyramid'.
    feat_f, faces_f, flags_f:
        Floating mesh arrays.
    feat_t, faces_t, flags_t:
        Target mesh arrays.
    """
    import meshmonk

    if mode == "rigid":
        meshmonk.rigid_register(
            floating_features=feat_f,
            target_features=feat_t,
            floating_faces=faces_f,
            target_faces=faces_t,
            floating_flags=flags_f,
            target_flags=flags_t,
        )
    elif mode == "nonrigid":
        meshmonk.nonrigid_register(
            floating_features=feat_f,
            target_features=feat_t,
            floating_faces=faces_f,
            target_faces=faces_t,
            floating_flags=flags_f,
            target_flags=flags_t,
        )
    elif mode == "pyramid":
        meshmonk.pyramid_register(
            floating_features=feat_f,
            target_features=feat_t,
            floating_faces=faces_f,
            target_faces=faces_t,
            floating_flags=flags_f,
            target_flags=flags_t,
        )
    else:
        raise ValueError(f"Unknown mode: {mode!r}")


def collect_stats(measured: list, overhead_us: float) -> dict:
    """Compute per-label statistics from a list of profiling_dump() results.

    Parameters
    ----------
    measured:
        List of profiling_dump() dicts, one per measured run.
        Each dict is ``{label: {"total_us": int, "count": int}}``.
    overhead_us:
        Empty-ScopedTimer overhead in microseconds (from profiling_calibrate).

    Returns
    -------
    dict
        ``{label: {median_total_us, median_count, corrected_total_us,
                   ms_per_call, share_pct, amdahl, fork_join_us, unreliable}}``
    """
    # Gather all labels across all measured runs
    all_labels = set()
    for run in measured:
        all_labels.update(run.keys())

    stats = {}
    for label in all_labels:
        total_us_vals = [run[label]["total_us"] for run in measured if label in run]
        count_vals = [run[label]["count"] for run in measured if label in run]

        median_total_us = float(np.median(total_us_vals))
        median_count = float(np.median(count_vals))

        # Corrected total: subtract timer overhead from each invocation
        corrected_total_us = median_total_us - median_count * overhead_us
        # Clamp to zero (overhead can exceed measured value for tiny functions)
        corrected_total_us = max(0.0, corrected_total_us)

        # ms per call (corrected)
        ms_per_call = (corrected_total_us / median_count / 1000.0) if median_count > 0 else 0.0

        stats[label] = {
            "median_total_us": median_total_us,
            "median_count": median_count,
            "corrected_total_us": corrected_total_us,
            "ms_per_call": ms_per_call,
            # share_pct and amdahl computed below (need total across all labels)
            "fork_join_us": 15.0 * median_count,
            # UNRELIABLE: ms_per_call < 5 * overhead_ms
            "unreliable": ms_per_call < 5.0 * (overhead_us / 1000.0),
        }

    # Compute share_pct and amdahl based on total corrected time
    total_corrected = sum(v["corrected_total_us"] for v in stats.values())
    for label, v in stats.items():
        if total_corrected > 0:
            share = v["corrected_total_us"] / total_corrected * 100.0
        else:
            share = 0.0
        v["share_pct"] = share
        # Amdahl ceiling: 1 / (1 - share/100), clamp denominator
        frac = share / 100.0
        if frac >= 1.0:
            v["amdahl"] = float("inf")
        elif frac <= 0.0:
            v["amdahl"] = 1.0
        else:
            v["amdahl"] = 1.0 / (1.0 - frac)

    return stats


def fit_scaling_exponent(
    stats_by_tier_mode: dict,
    label: str,
    mode: str,
) -> float:
    """Fit a log-log scaling exponent for a label across 1K/10K/100K tiers.

    Uses numpy.polyfit(log10(ns), log10(ys), 1) where ns are vertex counts
    and ys are corrected total times.

    Parameters
    ----------
    stats_by_tier_mode:
        Dict keyed by (tier, mode), values are ``{label: stats_dict}``.
    label:
        The label to fit.
    mode:
        The registration mode to use.

    Returns
    -------
    float
        Slope k from log-log fit, or float('nan') if insufficient data.
    """
    tier_n_map = {"1k": 1_000, "10k": 10_000, "100k": 100_000}
    tiers = ["1k", "10k", "100k"]

    xs = []
    ys = []
    for tier in tiers:
        key = (tier, mode)
        if key not in stats_by_tier_mode:
            continue
        label_stats = stats_by_tier_mode[key]
        if label not in label_stats:
            continue
        val = label_stats[label].get("corrected_total_us", 0.0)
        if val <= 0:
            continue
        xs.append(tier_n_map[tier])
        ys.append(val)

    if len(xs) < 2:
        return float("nan")

    log_xs = np.log10(xs)
    log_ys = np.log10(ys)
    coeffs = np.polyfit(log_xs, log_ys, 1)
    return float(coeffs[0])


def _complexity_bucket(k: float) -> str:
    """Assign a complexity bucket label given scaling exponent k."""
    if np.isnan(k):
        return "n/a (insufficient tiers)"
    candidates = [1.0, 1.16, 2.0]
    for c in candidates:
        if abs(k - c) < 0.2:
            return f"≈O(n^{c})"
    return f"other (k={k:.2f})"


def _pyramid_reconcile_check(stats: dict, mode: str) -> str:
    """Check that bare NonrigidRegistration::update ≈ sum of /layerN labels.

    Returns a warning string if divergence exceeds 5%, empty string otherwise.
    """
    if mode != "pyramid":
        return ""

    bare_key = "NonrigidRegistration::update"
    layer_keys = [k for k in stats if k.startswith("NonrigidRegistration::update/layer")]

    if not layer_keys or bare_key not in stats:
        return ""

    bare_total = stats[bare_key]["corrected_total_us"]
    layer_sum = sum(stats[k]["corrected_total_us"] for k in layer_keys)

    if layer_sum <= 0:
        return ""

    divergence = abs(bare_total - layer_sum) / layer_sum
    if divergence > 0.05:
        return (
            f"  > **SANITY WARNING**: bare `NonrigidRegistration::update` total "
            f"({bare_total:.0f} us) diverges from layer sum ({layer_sum:.0f} us) "
            f"by {divergence * 100:.1f}% (>5% tolerance).\n"
        )
    return ""


def generate_report(
    args,
    overhead_ns: int,
    stats_by_tier_mode: dict,
    scaling_exponents: dict,
) -> str:
    """Generate the full markdown profiling report.

    Parameters
    ----------
    args:
        Parsed argparse.Namespace from parse_args().
    overhead_ns:
        Empty-ScopedTimer overhead in nanoseconds (from profiling_calibrate).
    stats_by_tier_mode:
        Dict keyed by (tier, mode) -> {label -> stats_dict}.
    scaling_exponents:
        Dict keyed by (label, mode) -> float k (or NaN).

    Returns
    -------
    str
        Full markdown report text.
    """
    today = datetime.date.today().isoformat()

    # --- System info ---
    cpu_model = platform.processor() or platform.machine() or "unknown"
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("model name"):
                    cpu_model = line.split(":", 1)[1].strip()
                    break
    except Exception:
        pass

    try:
        import multiprocessing
        cpu_cores = multiprocessing.cpu_count()
    except Exception:
        cpu_cores = "unknown"

    os_info = f"{platform.system()} {platform.release()}"

    try:
        compiler = subprocess.check_output(["cc", "--version"]).decode().splitlines()[0]
    except Exception:
        compiler = "unknown"

    try:
        git_sha = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd="/workspace"
        ).strip().decode()
        git_short = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd="/workspace"
        ).strip().decode()
    except Exception:
        git_sha = "unknown"
        git_short = "unknown"

    # Mesh SHA-256 for each tier's Template + DemoFace
    mesh_sha = {}
    for name in [
        "Template_1K", "Template_10K", "Template_100K",
        "DemoFace_1K", "DemoFace_10K", "DemoFace_100K",
    ]:
        p = Path(f"/workspace/data/{name}.obj")
        if p.exists():
            mesh_sha[name] = sha256_file(p)
        else:
            mesh_sha[name] = "file not found"

    run_cmd = (
        f"python -m profiling.run_profile"
        f" --tiers {','.join(args.tiers)}"
        f" --modes {','.join(args.modes)}"
        f" --runs {args.runs}"
        f" --warmup {args.warmup}"
        f" --seed {args.seed}"
    )

    # -----------------------------------------------------------------------
    # Section 1: Methodology
    # -----------------------------------------------------------------------
    mesh_sha_rows = "\n".join(
        f"| {name}.obj | `{sha}` |" for name, sha in mesh_sha.items()
    )

    methodology = f"""# MeshMonk Hotspot Profile — {today} (`{git_short}`)

## Methodology

- **CPU**: {cpu_model} ({cpu_cores} logical processors)
- **OS**: {os_info}
- **Compiler**: {compiler}
- **CMAKE_BUILD_TYPE**: Release (assumed; use `cmake --build . --config Release`)
- **LTO**: unknown
- **Eigen**: vendored (vendor/eigen-3.4.0)
- **nanoflann**: vendored (vendor/nanoflann.hpp)
- **OpenMesh**: vendored (vendor/OpenMesh-11.0.0)
- **Repo SHA**: `{git_sha}`
- **OMP_NUM_THREADS**: {os.environ.get("OMP_NUM_THREADS", "not set")}
- **EIGEN_DONT_PARALLELIZE**: {os.environ.get("EIGEN_DONT_PARALLELIZE", "not set")}
- **Run command**: `{run_cmd}`
- **Seed**: {args.seed}
- **ScopedTimer overhead**: {overhead_ns} ns/scope (measured via profiling_calibrate(1_000_000))
- **Repetitions per (tier, mode)**: warmup={args.warmup}, runs={args.runs}, median selected

### Mesh SHA-256

| File | SHA-256 |
|---|---|
{mesh_sha_rows}

### Pyramid-mode deduplication note

In pyramid mode, `NonrigidRegistration::update` accumulates across all layers
(bare label, count = num_layers), AND each layer also emits a
`NonrigidRegistration::update/layerN` label. The bare label is included in
the table as a cross-check. Share percentages are computed from all labels
including the bare one — do not double-count by adding bare + layer shares.
A sanity check verifies bare_total ≈ sum(layerN_totals) within ±5%; violations
are noted inline.
"""

    # -----------------------------------------------------------------------
    # Section 2: Per-tier ranked tables
    # -----------------------------------------------------------------------
    n_verts_map = {"1k": "~1,000", "10k": "~10,000", "100k": "~100,000"}

    per_tier_sections = ["## Per-tier Hotspot Rankings\n"]

    for tier in args.tiers:
        for mode in args.modes:
            key = (tier, mode)
            if key not in stats_by_tier_mode:
                continue
            stats = stats_by_tier_mode[key]

            # Pyramid reconciliation check
            reconcile_warning = _pyramid_reconcile_check(stats, mode)

            # Sort by share_pct descending
            sorted_labels = sorted(
                stats.keys(), key=lambda lbl: stats[lbl]["share_pct"], reverse=True
            )

            rows = []
            for label in sorted_labels:
                v = stats[label]
                notes = "* UNRELIABLE" if v["unreliable"] else ""
                rows.append(
                    f"| `{label}` "
                    f"| {v['share_pct']:.1f} "
                    f"| {v['ms_per_call']:.4f} "
                    f"| {int(v['median_count'])} "
                    f"| {v['amdahl']:.1f}x "
                    f"| {v['fork_join_us']:.0f} "
                    f"| {notes} |"
                )

            table_rows = "\n".join(rows) if rows else "| (no data) | — | — | — | — | — | |"

            per_tier_sections.append(
                f"\n### Tier: {tier.upper()} ({n_verts_map.get(tier, '?')} vertices) — {mode}\n\n"
                f"| Label | Share (%) | ms/call | Invocations | Amdahl ceiling"
                f" | Fork-join cost (us) | Notes |\n"
                f"|---|---|---|---|---|---|---|\n"
                f"{table_rows}\n\n"
                f"{reconcile_warning}"
                f"*UNRELIABLE = ms/call < 5× timer overhead"
                f" ({overhead_ns} ns/scope = {overhead_ns / 1e6:.4f} ms)."
                f" Excluded from recommendations.*\n"
            )

    # -----------------------------------------------------------------------
    # Section 3: Scaling exponents
    # -----------------------------------------------------------------------
    scaling_rows = []
    seen_label_mode = set()
    for tier in args.tiers:
        for mode in args.modes:
            key = (tier, mode)
            if key not in stats_by_tier_mode:
                continue
            for label in stats_by_tier_mode[key]:
                lm_key = (label, mode)
                if lm_key in seen_label_mode:
                    continue
                seen_label_mode.add(lm_key)
                k = scaling_exponents.get(lm_key, float("nan"))
                bucket = _complexity_bucket(k)
                k_str = f"{k:.2f}" if not np.isnan(k) else "NaN"
                note = "insufficient tiers" if np.isnan(k) else ""
                scaling_rows.append(
                    f"| `{label}` | {mode} | {k_str} | {bucket} | {note} |"
                )

    scaling_table = "\n".join(scaling_rows) if scaling_rows else "| (no data) | — | — | — | — |"

    scaling_section = f"""## Scaling Exponents (log-log fit across 1K / 10K / 100K)

| Label | Mode | k (exponent) | Complexity bucket | Note |
|---|---|---|---|---|
{scaling_table}

> **Caveat:** Three-point fits (1K/10K/100K) have zero residual degrees of freedom.
> The buckets k≈1.0 / k≈1.16 / k≈2.0 cannot be reliably discriminated without a
> fourth tier or per-tier variance. Treat exponents as indicative only.
"""

    # -----------------------------------------------------------------------
    # Section 4: Recommendations
    # -----------------------------------------------------------------------
    candidates = []
    below_threshold = []
    excluded_unreliable = []

    for tier in args.tiers:
        for mode in args.modes:
            key = (tier, mode)
            if key not in stats_by_tier_mode:
                continue
            stats = stats_by_tier_mode[key]
            for label, v in stats.items():
                if v["unreliable"]:
                    excluded_unreliable.append((label, mode, tier, v))
                elif v["share_pct"] > 5.0 and v["ms_per_call"] > 0.03:
                    candidates.append((label, mode, tier, v))
                else:
                    below_threshold.append((label, mode, tier, v))

    # Deduplicate candidates by (label, mode) keeping highest share
    best_candidates = {}
    for label, mode, tier, v in candidates:
        lm = (label, mode)
        if lm not in best_candidates or v["share_pct"] > best_candidates[lm][3]["share_pct"]:
            best_candidates[lm] = (label, mode, tier, v)
    candidates = sorted(best_candidates.values(), key=lambda x: x[3]["share_pct"], reverse=True)

    candidate_text_parts = []
    for label, mode, tier, v in candidates:
        k = scaling_exponents.get((label, mode), float("nan"))
        k_str = f"{k:.2f}" if not np.isnan(k) else "NaN"
        # Tool suggestion heuristic
        fj_ratio = v["fork_join_us"] / (v["corrected_total_us"] + 1e-9) * 100.0
        if v["ms_per_call"] > 5.0 * (15e-3):  # fork-join ~15us threshold
            tool_suggest = "OpenMP (ms/call >> fork-join cost)"
        else:
            tool_suggest = "algorithmic / Eigen-level (fork-join overhead too high for OpenMP)"
        fj_warn = ""
        if fj_ratio > 10.0:
            fj_warn = (
                f"\n  > **NOTE**: Fork-join cost is {fj_ratio:.1f}% of label wall time —"
                f" OpenMP overhead likely exceeds speedup."
            )
        candidate_text_parts.append(
            f"- **`{label}`** ({mode}, {tier.upper()} tier):"
            f" share={v['share_pct']:.1f}%,"
            f" Amdahl ceiling = 1/(1−{v['share_pct'] / 100:.3f}) = {v['amdahl']:.1f}×\n"
            f"  - Scaling exponent k: {k_str}\n"
            f"  - Suggested tool: {tool_suggest}\n"
            f"  - Fork-join cost: {v['fork_join_us']:.0f} µs"
            f" ({fj_ratio:.1f}% of label wall time){fj_warn}"
        )

    candidates_md = "\n\n".join(candidate_text_parts) if candidate_text_parts else "(none)"

    # Below threshold summary
    below_labels = set(label for label, _, _, _ in below_threshold)
    below_md = (
        ", ".join(f"`{lbl}`" for lbl in sorted(below_labels))
        if below_labels
        else "(none)"
    )

    # Unreliable summary
    unreliable_rows = []
    seen_ur = set()
    for label, mode, tier, v in excluded_unreliable:
        if (label, mode) in seen_ur:
            continue
        seen_ur.add((label, mode))
        unreliable_rows.append(
            f"| `{label}` | {mode} | {tier.upper()} | {v['ms_per_call']:.4f} |"
        )
    unreliable_table = (
        "| Label | Mode | Tier | ms/call |\n|---|---|---|---|\n"
        + "\n".join(unreliable_rows)
        if unreliable_rows
        else "(none)"
    )

    recommendation_section = f"""## Recommendations

### Candidates for optimization

Labels meeting all of: share > 5%, ms/call > 0.03 ms, NOT UNRELIABLE.

{candidates_md}

### Below granularity threshold

Labels with share ≤ 5% or ms/call ≤ 0.03 ms. OpenMP not viable due to
fork-join overhead exceeding potential speedup (reference: beads 9f5, bdt).

{below_md}

### Excluded (UNRELIABLE)

Labels with ms/call < 5× timer overhead ({overhead_ns} ns/scope = {overhead_ns / 1e6:.4f} ms).
Corrected values unreliable — cannot meaningfully distinguish measurement noise from actual cost.

{unreliable_table}
"""

    # -----------------------------------------------------------------------
    # Assemble full report
    # -----------------------------------------------------------------------
    return (
        methodology
        + "\n"
        + "\n".join(per_tier_sections)
        + "\n"
        + scaling_section
        + "\n"
        + recommendation_section
    )


def main() -> None:
    """Main entry point for the profiling driver."""
    # Step 0: Pin env vars BEFORE importing meshmonk — C++ reads them at load time
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["EIGEN_DONT_PARALLELIZE"] = "1"
    import meshmonk  # deferred import — env vars must be set first

    args = parse_args()

    # Step 1: Preflight assertion
    if not meshmonk.profiling_enabled():
        print(
            "ERROR: meshmonk was not built with -DMESHMONK_PROFILING=ON "
            "-- the report would be all zeros",
            file=sys.stderr,
        )
        sys.exit(2)

    # Step 2: Calibration
    N_CAL = 1_000_000
    overhead_ns = meshmonk.profiling_calibrate(N_CAL)
    overhead_us = overhead_ns / 1000.0
    print(f"ScopedTimer overhead: {overhead_ns} ns/scope ({overhead_us:.3f} us/scope)")

    # Step 3 + 4 + 5: Mesh loading and measurement loop
    stats_by_tier_mode = {}

    for tier in args.tiers:
        print(f"\nLoading meshes for tier {tier.upper()} ...")
        feat_f, faces_f, flags_f, feat_t, faces_t, flags_t = load_mesh_tier(tier)
        n_verts = feat_f.shape[0]
        print(f"  Floating: {n_verts} vertices, Target: {feat_t.shape[0]} vertices")

        for mode in args.modes:
            print(f"  Running {mode} at {tier.upper()} ({args.warmup} warmup + {args.runs} runs) ...")
            all_data = []
            for rep in range(args.warmup + args.runs):
                meshmonk.profiling_reset()
                run_registration(mode, feat_f, faces_f, flags_f, feat_t, faces_t, flags_t)
                data = meshmonk.profiling_dump()
                all_data.append(data)
                if rep < args.warmup:
                    print(f"    rep {rep + 1}/{args.warmup + args.runs} (warmup) — {len(data)} labels")
                else:
                    print(f"    rep {rep + 1}/{args.warmup + args.runs} — {len(data)} labels")
            measured = all_data[args.warmup:]  # discard warmup runs
            stats_by_tier_mode[(tier, mode)] = collect_stats(measured, overhead_us)
            n_labels = len(stats_by_tier_mode[(tier, mode)])
            print(f"  {mode} at {tier.upper()}: {n_labels} labels collected")

    # Step 6: Compute scaling exponents
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

    # Step 7: Generate report
    print("\nGenerating report ...")
    report_md = generate_report(args, overhead_ns, stats_by_tier_mode, scaling_exponents)

    # Step 8: Write report to file
    if args.out:
        out_path = Path(args.out)
    else:
        try:
            git_short = subprocess.check_output(
                ["git", "rev-parse", "--short", "HEAD"], cwd="/workspace"
            ).strip().decode()
        except Exception:
            git_short = "unknown"
        today = datetime.date.today().strftime("%Y%m%d")
        out_dir = Path("/workspace/docs/perf")
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"hotspot-profile-{today}-{git_short}.md"

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(report_md, encoding="utf-8")
    print(f"\nReport written to: {out_path}")

    # Step 9: Sanity checks
    print("\nRunning sanity checks ...")
    report_text = report_md

    checks = [
        ("NeighbourFinder::update", "NeighbourFinder::update missing from report"),
        ("ViscoElasticTransformer", "No ViscoElasticTransformer label in report"),
        ("NonrigidRegistration::update", "NonrigidRegistration::update missing"),
        ("RigidRegistration::update", "RigidRegistration::update missing"),
        ("OMP_NUM_THREADS", "OMP_NUM_THREADS not in Methodology"),
        ("ns/scope", "Calibration overhead not in Methodology"),
    ]

    # Only check pyramid labels if pyramid was run
    if "pyramid" in args.modes:
        checks.append(
            ("PyramidNonrigidRegistration::update", "Pyramid static label missing")
        )
        checks.append(("/layer", "No /layer labels in pyramid mode"))

    all_passed = True
    for needle, msg in checks:
        if needle not in report_text:
            print(f"  FAIL: {msg}")
            all_passed = False
        else:
            print(f"  PASS: '{needle}' found")

    if all_passed:
        print("\nAll sanity checks passed — report is ready.")
    else:
        print("\nWARNING: Some sanity checks failed — review the report before committing.")
        sys.exit(1)


if __name__ == "__main__":
    main()
