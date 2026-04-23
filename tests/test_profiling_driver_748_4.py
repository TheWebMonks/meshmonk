"""Tests for the profiling driver script (bead meshmonk-modernization-748.4).

Tests are structured to cover:
- parse_args: argument parsing and defaults
- sha256_file: deterministic SHA-256 hash
- load_mesh_tier: loads all 6 OBJ files correctly
- run_registration: dispatches to correct meshmonk function per mode
- collect_stats: computes corrected stats, unreliable flags
- fit_scaling_exponent: log-log regression, NaN on insufficient data
- generate_report: produces markdown with all 4 required sections
- main: step-0 guard (profiling_enabled=False -> sys.exit(2))

All tests are designed to work when MESHMONK_PROFILING is ON.
Tests that require profiling_enabled() are skipped otherwise.
"""

import sys
from pathlib import Path
from unittest.mock import patch, MagicMock
import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Import the module under test
# ---------------------------------------------------------------------------

sys.path.insert(0, str(Path("/workspace")))
import profiling.run_profile as rp


# ---------------------------------------------------------------------------
# parse_args
# ---------------------------------------------------------------------------

class TestParseArgs:
    def test_defaults(self):
        args = rp.parse_args(["--tiers", "1k", "--modes", "rigid"])
        assert args.tiers == ["1k"]
        assert args.modes == ["rigid"]
        assert args.runs == 5
        assert args.warmup == 1
        assert args.seed == 42

    def test_multiple_tiers_modes(self):
        args = rp.parse_args(["--tiers", "1k,10k,100k", "--modes", "rigid,nonrigid,pyramid"])
        assert args.tiers == ["1k", "10k", "100k"]
        assert args.modes == ["rigid", "nonrigid", "pyramid"]

    def test_custom_runs_warmup(self):
        args = rp.parse_args(["--tiers", "1k", "--modes", "rigid", "--runs", "3", "--warmup", "2"])
        assert args.runs == 3
        assert args.warmup == 2

    def test_out_argument(self):
        args = rp.parse_args(["--tiers", "1k", "--modes", "rigid", "--out", "foo.md"])
        assert args.out == "foo.md"

    def test_out_defaults_to_none(self):
        args = rp.parse_args(["--tiers", "1k", "--modes", "rigid"])
        assert args.out is None


# ---------------------------------------------------------------------------
# sha256_file
# ---------------------------------------------------------------------------

class TestSha256File:
    def test_consistent_hash(self, tmp_path):
        f = tmp_path / "test.txt"
        f.write_bytes(b"hello world")
        h1 = rp.sha256_file(f)
        h2 = rp.sha256_file(f)
        assert h1 == h2

    def test_different_content_different_hash(self, tmp_path):
        f1 = tmp_path / "a.txt"
        f2 = tmp_path / "b.txt"
        f1.write_bytes(b"content A")
        f2.write_bytes(b"content B")
        assert rp.sha256_file(f1) != rp.sha256_file(f2)

    def test_hash_is_hex_string(self, tmp_path):
        f = tmp_path / "test.txt"
        f.write_bytes(b"data")
        h = rp.sha256_file(f)
        assert isinstance(h, str)
        assert len(h) == 64
        int(h, 16)  # valid hex


# ---------------------------------------------------------------------------
# load_mesh_tier
# ---------------------------------------------------------------------------

class TestLoadMeshTier:
    """Requires the 6 OBJ files in /workspace/data/"""

    DATA = Path("/workspace/data")

    def test_1k_returns_correct_shapes(self):
        if not (self.DATA / "Template_1K.obj").exists():
            pytest.skip("1K mesh files not found")
        feat_f, ff, flags_f, feat_t, ft, flags_t = rp.load_mesh_tier("1k")
        # features should be (N,6) float32
        assert feat_f.ndim == 2
        assert feat_f.shape[1] == 6
        assert feat_f.dtype == np.float32
        # faces should be (M,3) int32
        assert ff.ndim == 2
        assert ff.shape[1] == 3
        assert ff.dtype == np.int32
        # flags should be (N,) float32
        assert flags_f.ndim == 1
        assert flags_f.dtype == np.float32
        assert len(flags_f) == feat_f.shape[0]

    def test_all_tiers_loadable(self):
        for tier in ["1k", "10k", "100k"]:
            tier_files = {
                "1k": ("Template_1K.obj", "DemoFace_1K.obj"),
                "10k": ("Template_10K.obj", "DemoFace_10K.obj"),
                "100k": ("Template_100K.obj", "DemoFace_100K.obj"),
            }
            t_file, d_file = tier_files[tier]
            if not (self.DATA / t_file).exists():
                pytest.skip(f"{t_file} not found")
            result = rp.load_mesh_tier(tier)
            assert len(result) == 6

    def test_tier_case_insensitive(self):
        if not (self.DATA / "Template_1K.obj").exists():
            pytest.skip("1K mesh files not found")
        result_lower = rp.load_mesh_tier("1k")
        result_upper = rp.load_mesh_tier("1K")
        assert result_lower[0].shape == result_upper[0].shape

    def test_flags_all_ones(self):
        if not (self.DATA / "Template_1K.obj").exists():
            pytest.skip("1K mesh files not found")
        feat_f, ff, flags_f, feat_t, ft, flags_t = rp.load_mesh_tier("1k")
        np.testing.assert_array_equal(flags_f, np.ones(len(flags_f)))
        np.testing.assert_array_equal(flags_t, np.ones(len(flags_t)))


# ---------------------------------------------------------------------------
# run_registration (mocked meshmonk)
# ---------------------------------------------------------------------------

class TestRunRegistration:
    """Tests run_registration dispatches correctly. Mocks meshmonk."""

    def _make_arrays(self, n=10):
        faces = np.zeros((n // 3, 3), dtype=np.int32)
        feat = np.random.rand(n, 6).astype(np.float32)
        flags = np.ones(n, dtype=np.float32)
        return feat, faces, flags

    def test_rigid_mode_calls_rigid_register(self):
        feat_f, ff, flags_f = self._make_arrays()
        feat_t, ft, flags_t = self._make_arrays()
        mock_mm = MagicMock()
        with patch.dict("sys.modules", {"meshmonk": mock_mm}):
            rp.run_registration("rigid", feat_f, ff, flags_f, feat_t, ft, flags_t)
        mock_mm.rigid_register.assert_called_once()

    def test_nonrigid_mode_calls_nonrigid_register(self):
        feat_f, ff, flags_f = self._make_arrays()
        feat_t, ft, flags_t = self._make_arrays()
        mock_mm = MagicMock()
        with patch.dict("sys.modules", {"meshmonk": mock_mm}):
            rp.run_registration("nonrigid", feat_f, ff, flags_f, feat_t, ft, flags_t)
        mock_mm.nonrigid_register.assert_called_once()

    def test_pyramid_mode_calls_pyramid_register(self):
        feat_f, ff, flags_f = self._make_arrays()
        feat_t, ft, flags_t = self._make_arrays()
        mock_mm = MagicMock()
        with patch.dict("sys.modules", {"meshmonk": mock_mm}):
            rp.run_registration("pyramid", feat_f, ff, flags_f, feat_t, ft, flags_t)
        mock_mm.pyramid_register.assert_called_once()

    def test_unknown_mode_raises(self):
        feat_f, ff, flags_f = self._make_arrays()
        feat_t, ft, flags_t = self._make_arrays()
        with pytest.raises(ValueError, match="Unknown mode"):
            rp.run_registration("foobar", feat_f, ff, flags_f, feat_t, ft, flags_t)

    def test_rigid_passes_correct_kwargs(self):
        feat_f, ff, flags_f = self._make_arrays()
        feat_t, ft, flags_t = self._make_arrays()
        mock_mm = MagicMock()
        with patch.dict("sys.modules", {"meshmonk": mock_mm}):
            rp.run_registration("rigid", feat_f, ff, flags_f, feat_t, ft, flags_t)
        call_kwargs = mock_mm.rigid_register.call_args.kwargs
        assert "floating_features" in call_kwargs
        assert "target_features" in call_kwargs
        assert "floating_faces" in call_kwargs
        assert "target_faces" in call_kwargs
        assert "floating_flags" in call_kwargs
        assert "target_flags" in call_kwargs


# ---------------------------------------------------------------------------
# collect_stats
# ---------------------------------------------------------------------------

class TestCollectStats:
    """Tests for the stats collection and correction logic."""

    def _make_measured(self, label, total_us_values, count_values):
        """Build a list of profiling_dump() dicts with consistent shape."""
        return [
            {label: {"total_us": t, "count": c}}
            for t, c in zip(total_us_values, count_values)
        ]

    def test_median_total_us_computed(self):
        measured = self._make_measured(
            "TestLabel::fn",
            total_us_values=[100, 200, 150, 180, 160],
            count_values=[10, 10, 10, 10, 10],
        )
        stats = rp.collect_stats(measured, overhead_us=0.0)
        assert "TestLabel::fn" in stats
        assert stats["TestLabel::fn"]["median_total_us"] == pytest.approx(160.0)

    def test_corrected_total_subtracts_overhead(self):
        # 5 runs, count=10, overhead_us=1.0
        # corrected = median_total - median_count * overhead_us
        measured = self._make_measured(
            "TestLabel::fn",
            total_us_values=[200, 200, 200, 200, 200],
            count_values=[10, 10, 10, 10, 10],
        )
        stats = rp.collect_stats(measured, overhead_us=1.0)
        # corrected = 200 - 10 * 1.0 = 190
        assert stats["TestLabel::fn"]["corrected_total_us"] == pytest.approx(190.0)

    def test_ms_per_call_computed(self):
        # corrected = 200 - 10*0 = 200 us total, count=10 -> 20 us/call -> 0.020 ms/call
        measured = self._make_measured(
            "TestLabel::fn",
            total_us_values=[200, 200, 200, 200, 200],
            count_values=[10, 10, 10, 10, 10],
        )
        stats = rp.collect_stats(measured, overhead_us=0.0)
        assert stats["TestLabel::fn"]["ms_per_call"] == pytest.approx(200.0 / 10.0 / 1000.0)

    def test_share_pct_sums_to_100_with_two_labels(self):
        measured = [
            {"A": {"total_us": 300, "count": 1}, "B": {"total_us": 700, "count": 1}}
            for _ in range(5)
        ]
        stats = rp.collect_stats(measured, overhead_us=0.0)
        total_share = stats["A"]["share_pct"] + stats["B"]["share_pct"]
        assert total_share == pytest.approx(100.0, abs=0.01)

    def test_unreliable_flag_when_ms_per_call_lt_5x_overhead(self):
        # overhead_us=1.0 -> overhead_ms=0.001, 5x = 0.005 ms
        # ms_per_call = 100us / 1 / 1000 = 0.1 ms -> NOT unreliable
        measured = self._make_measured(
            "FastLabel",
            total_us_values=[100, 100, 100, 100, 100],
            count_values=[1, 1, 1, 1, 1],
        )
        stats = rp.collect_stats(measured, overhead_us=0.0)
        assert not stats["FastLabel"]["unreliable"]

    def test_unreliable_flag_set_when_ms_per_call_lt_5x_overhead(self):
        # overhead_us = 10.0 -> overhead_ms = 0.01 ms, 5x = 0.05 ms
        # corrected = 20 us/call -> 0.02 ms/call < 5x overhead -> UNRELIABLE
        measured = self._make_measured(
            "TinyLabel",
            total_us_values=[20, 20, 20, 20, 20],
            count_values=[1, 1, 1, 1, 1],
        )
        stats = rp.collect_stats(measured, overhead_us=10.0)
        assert stats["TinyLabel"]["unreliable"]

    def test_amdahl_computed(self):
        # share 50% -> amdahl = 1/(1-0.5) = 2.0x
        measured = [
            {"A": {"total_us": 500, "count": 1}, "B": {"total_us": 500, "count": 1}}
            for _ in range(5)
        ]
        stats = rp.collect_stats(measured, overhead_us=0.0)
        assert stats["A"]["amdahl"] == pytest.approx(2.0, rel=0.01)

    def test_fork_join_cost_computed(self):
        # fork_join_us = 15 * count
        measured = self._make_measured(
            "Func",
            total_us_values=[1000, 1000, 1000, 1000, 1000],
            count_values=[20, 20, 20, 20, 20],
        )
        stats = rp.collect_stats(measured, overhead_us=0.0)
        assert stats["Func"]["fork_join_us"] == pytest.approx(15 * 20)

    def test_labels_from_all_runs_collected(self):
        # label present in some but not all runs
        measured = [
            {"A": {"total_us": 100, "count": 1}},
            {"A": {"total_us": 110, "count": 1}, "B": {"total_us": 50, "count": 1}},
            {"A": {"total_us": 105, "count": 1}, "B": {"total_us": 55, "count": 1}},
            {"A": {"total_us": 100, "count": 1}},
            {"A": {"total_us": 102, "count": 1}, "B": {"total_us": 52, "count": 1}},
        ]
        stats = rp.collect_stats(measured, overhead_us=0.0)
        assert "A" in stats
        assert "B" in stats

    def test_pyramid_bare_label_reconciliation(self):
        """In pyramid mode, bare NonrigidRegistration::update should match sum of layers.
        collect_stats stores the raw labels - reconciliation is done in generate_report."""
        # Simulate data like what pyramid produces
        layer_times = [24521, 24954, 28931]
        bare_total = sum(layer_times)  # Perfect match scenario
        measured = [
            {
                "NonrigidRegistration::update": {"total_us": bare_total, "count": 3},
                "NonrigidRegistration::update/layer0": {"total_us": layer_times[0], "count": 1},
                "NonrigidRegistration::update/layer1": {"total_us": layer_times[1], "count": 1},
                "NonrigidRegistration::update/layer2": {"total_us": layer_times[2], "count": 1},
                "PyramidNonrigidRegistration::update": {"total_us": bare_total + 5000, "count": 1},
            }
            for _ in range(5)
        ]
        stats = rp.collect_stats(measured, overhead_us=0.0)
        # All labels present
        assert "NonrigidRegistration::update" in stats
        assert "NonrigidRegistration::update/layer0" in stats
        assert "PyramidNonrigidRegistration::update" in stats


# ---------------------------------------------------------------------------
# fit_scaling_exponent
# ---------------------------------------------------------------------------

class TestFitScalingExponent:
    """Tests for the log-log regression across tiers."""

    def _make_stats_by_tier_mode(self, label, mode, values):
        """Build stats_by_tier_mode dict with the given total_us values per tier."""
        tiers = ["1k", "10k", "100k"]
        result = {}
        for tier, val in zip(tiers, values):
            result[(tier, mode)] = {
                label: {
                    "corrected_total_us": val,
                    "median_count": 1,
                    "ms_per_call": val / 1000.0,
                }
            }
        return result

    def test_linear_scaling_gives_k1(self):
        # If cost scales linearly with n, k should be ~1.0
        # 1k->10x->100x means y1, y2=10*y1, y3=100*y1
        stats = self._make_stats_by_tier_mode("Func", "rigid", [100, 1000, 10000])
        k = rp.fit_scaling_exponent(stats, "Func", "rigid")
        assert k == pytest.approx(1.0, abs=0.05)

    def test_quadratic_scaling_gives_k2(self):
        # Quadratic: y values grow as n^2
        # n=1k->1, n=10k->100, n=100k->10000
        stats = self._make_stats_by_tier_mode("Func", "rigid", [1, 100, 10000])
        k = rp.fit_scaling_exponent(stats, "Func", "rigid")
        assert k == pytest.approx(2.0, abs=0.05)

    def test_missing_label_returns_nan(self):
        # Label not present in any tier -> NaN
        stats = self._make_stats_by_tier_mode("OtherLabel", "rigid", [100, 1000, 10000])
        k = rp.fit_scaling_exponent(stats, "MissingLabel", "rigid")
        assert np.isnan(k)

    def test_label_in_only_one_tier_returns_nan(self):
        tiers = ["1k", "10k", "100k"]
        stats = {}
        for i, tier in enumerate(tiers):
            if i == 0:
                stats[(tier, "rigid")] = {
                    "Func": {"corrected_total_us": 100, "median_count": 1}
                }
            else:
                stats[(tier, "rigid")] = {}
        k = rp.fit_scaling_exponent(stats, "Func", "rigid")
        assert np.isnan(k)

    def test_nlogn_scaling_gives_k_approx_116(self):
        # n*log(n) scaling: approximate fit
        import math
        ns = [1000, 10000, 100000]
        ys = [n * math.log(n) for n in ns]
        stats = self._make_stats_by_tier_mode("Func", "rigid", ys)
        k = rp.fit_scaling_exponent(stats, "Func", "rigid")
        # Should be between 1.0 and 1.25 for n*log(n)
        assert 1.0 <= k <= 1.3


# ---------------------------------------------------------------------------
# generate_report
# ---------------------------------------------------------------------------

class TestGenerateReport:
    """Tests for the markdown report generation."""

    def _make_minimal_args(self):
        args = rp.parse_args(["--tiers", "1k", "--modes", "rigid"])
        return args

    def _make_stats_by_tier_mode(self):
        return {
            ("1k", "rigid"): {
                "RigidRegistration::update": {
                    "median_total_us": 50000.0,
                    "median_count": 1.0,
                    "corrected_total_us": 49900.0,
                    "ms_per_call": 49.9,
                    "share_pct": 74.5,
                    "amdahl": 3.9,
                    "fork_join_us": 15.0,
                    "unreliable": False,
                },
                "NeighbourFinder::update": {
                    "median_total_us": 17000.0,
                    "median_count": 10.0,
                    "corrected_total_us": 16900.0,
                    "ms_per_call": 1.69,
                    "share_pct": 25.2,
                    "amdahl": 1.34,
                    "fork_join_us": 150.0,
                    "unreliable": False,
                },
            }
        }

    def _make_scaling_exponents(self):
        return {("RigidRegistration::update", "rigid"): 1.05}

    def test_report_contains_methodology_section(self):
        report = rp.generate_report(
            self._make_minimal_args(),
            overhead_ns=100,
            stats_by_tier_mode=self._make_stats_by_tier_mode(),
            scaling_exponents=self._make_scaling_exponents(),
        )
        assert "## Methodology" in report

    def test_report_contains_per_tier_section(self):
        report = rp.generate_report(
            self._make_minimal_args(),
            overhead_ns=100,
            stats_by_tier_mode=self._make_stats_by_tier_mode(),
            scaling_exponents=self._make_scaling_exponents(),
        )
        assert "## Per-tier" in report

    def test_report_contains_scaling_section(self):
        report = rp.generate_report(
            self._make_minimal_args(),
            overhead_ns=100,
            stats_by_tier_mode=self._make_stats_by_tier_mode(),
            scaling_exponents=self._make_scaling_exponents(),
        )
        assert "## Scaling" in report

    def test_report_contains_recommendation_section(self):
        report = rp.generate_report(
            self._make_minimal_args(),
            overhead_ns=100,
            stats_by_tier_mode=self._make_stats_by_tier_mode(),
            scaling_exponents=self._make_scaling_exponents(),
        )
        assert "## Recommendation" in report

    def test_report_has_four_sections(self):
        report = rp.generate_report(
            self._make_minimal_args(),
            overhead_ns=100,
            stats_by_tier_mode=self._make_stats_by_tier_mode(),
            scaling_exponents=self._make_scaling_exponents(),
        )
        # Count top-level sections
        import re
        sections = re.findall(r"^## ", report, re.MULTILINE)
        assert len(sections) >= 4

    def test_report_contains_omp_threads(self):
        import os
        os.environ["OMP_NUM_THREADS"] = "1"
        report = rp.generate_report(
            self._make_minimal_args(),
            overhead_ns=100,
            stats_by_tier_mode=self._make_stats_by_tier_mode(),
            scaling_exponents=self._make_scaling_exponents(),
        )
        assert "OMP_NUM_THREADS" in report

    def test_report_contains_calibration_overhead(self):
        report = rp.generate_report(
            self._make_minimal_args(),
            overhead_ns=150,
            stats_by_tier_mode=self._make_stats_by_tier_mode(),
            scaling_exponents=self._make_scaling_exponents(),
        )
        assert "ns/scope" in report

    def test_report_contains_label_name(self):
        report = rp.generate_report(
            self._make_minimal_args(),
            overhead_ns=100,
            stats_by_tier_mode=self._make_stats_by_tier_mode(),
            scaling_exponents=self._make_scaling_exponents(),
        )
        assert "RigidRegistration::update" in report

    def test_report_unreliable_label_flagged(self):
        stats = {
            ("1k", "rigid"): {
                "TinyFunc": {
                    "median_total_us": 5.0,
                    "median_count": 1.0,
                    "corrected_total_us": 5.0,
                    "ms_per_call": 0.005,
                    "share_pct": 100.0,
                    "amdahl": 999.0,
                    "fork_join_us": 15.0,
                    "unreliable": True,
                }
            }
        }
        report = rp.generate_report(
            self._make_minimal_args(),
            overhead_ns=100,
            stats_by_tier_mode=stats,
            scaling_exponents={},
        )
        assert "UNRELIABLE" in report

    def test_pyramid_bare_label_reconciliation_note(self):
        """Pyramid mode report should note the dedup logic."""
        args = rp.parse_args(["--tiers", "1k", "--modes", "pyramid"])
        bare_total = 71366  # differs from layer sum (78406) by ~8.5%
        stats = {
            ("1k", "pyramid"): {
                "PyramidNonrigidRegistration::update": {
                    "median_total_us": 85000.0, "median_count": 1.0,
                    "corrected_total_us": 85000.0, "ms_per_call": 85.0,
                    "share_pct": 100.0, "amdahl": 999.0, "fork_join_us": 15.0,
                    "unreliable": False,
                },
                "NonrigidRegistration::update": {
                    "median_total_us": float(bare_total), "median_count": 3.0,
                    "corrected_total_us": float(bare_total), "ms_per_call": 23.79,
                    "share_pct": 84.0, "amdahl": 6.25, "fork_join_us": 45.0,
                    "unreliable": False,
                },
                "NonrigidRegistration::update/layer0": {
                    "median_total_us": 24521.0, "median_count": 1.0,
                    "corrected_total_us": 24521.0, "ms_per_call": 24.521,
                    "share_pct": 28.8, "amdahl": 1.4, "fork_join_us": 15.0,
                    "unreliable": False,
                },
                "NonrigidRegistration::update/layer1": {
                    "median_total_us": 24954.0, "median_count": 1.0,
                    "corrected_total_us": 24954.0, "ms_per_call": 24.954,
                    "share_pct": 29.3, "amdahl": 1.4, "fork_join_us": 15.0,
                    "unreliable": False,
                },
                "NonrigidRegistration::update/layer2": {
                    "median_total_us": 28931.0, "median_count": 1.0,
                    "corrected_total_us": 28931.0, "ms_per_call": 28.931,
                    "share_pct": 34.0, "amdahl": 1.5, "fork_join_us": 15.0,
                    "unreliable": False,
                },
            }
        }
        report = rp.generate_report(
            args, overhead_ns=100, stats_by_tier_mode=stats, scaling_exponents={}
        )
        # Should contain layer labels
        assert "/layer" in report
        # Should contain the pyramid outer label
        assert "PyramidNonrigidRegistration::update" in report


# ---------------------------------------------------------------------------
# main(): step-0 guard
# ---------------------------------------------------------------------------

class TestMainStepZeroGuard:
    """main() must exit with code 2 when profiling_enabled() returns False."""

    def test_exits_with_code_2_when_profiling_disabled(self):
        """Call rp.main() with mocked args and meshmonk; assert exit code 2."""
        mock_mm = MagicMock()
        mock_mm.profiling_enabled.return_value = False

        mock_args = MagicMock()
        mock_args.tiers = ["1k"]
        mock_args.modes = ["rigid"]
        mock_args.runs = 1
        mock_args.warmup = 0
        mock_args.seed = 42
        mock_args.out = None

        with patch.dict("sys.modules", {"meshmonk": mock_mm}):
            with patch.object(rp, "parse_args", return_value=mock_args):
                with pytest.raises(SystemExit) as exc_info:
                    rp.main()

        assert exc_info.value.code == 2


# ---------------------------------------------------------------------------
# Integration: import sanity
# ---------------------------------------------------------------------------

class TestModuleStructure:
    """Verify module-level structure."""

    def test_all_required_functions_exist(self):
        for fn in [
            "parse_args",
            "sha256_file",
            "load_mesh_tier",
            "run_registration",
            "collect_stats",
            "fit_scaling_exponent",
            "generate_report",
            "main",
        ]:
            assert hasattr(rp, fn), f"Missing function: {fn}"

    def test_parse_args_is_callable(self):
        assert callable(rp.parse_args)

    def test_generate_report_is_callable(self):
        assert callable(rp.generate_report)
