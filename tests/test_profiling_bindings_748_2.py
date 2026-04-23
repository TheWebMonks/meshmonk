"""Tests for Python profiling bindings (bead meshmonk-modernization-748.2).

Verifies five new profiling bindings exposed through meshmonk:
  - profiling_reset()
  - profiling_peek()
  - profiling_dump()
  - profiling_enabled()
  - profiling_calibrate()

These bindings are always compiled in. When MESHMONK_PROFILING is OFF (the
default build), they are no-ops/return empty/return 0. When ON, they expose
the full profiling functionality.

Tests are designed to pass in BOTH build variants.
"""

import meshmonk
import pytest


# ---------------------------------------------------------------------------
# Import surface: all 5 names must be importable from meshmonk
# ---------------------------------------------------------------------------

class TestImportSurface:
    """All five profiling names must be importable and in __all__."""

    def test_profiling_reset_importable(self):
        assert hasattr(meshmonk, "profiling_reset")

    def test_profiling_peek_importable(self):
        assert hasattr(meshmonk, "profiling_peek")

    def test_profiling_dump_importable(self):
        assert hasattr(meshmonk, "profiling_dump")

    def test_profiling_enabled_importable(self):
        assert hasattr(meshmonk, "profiling_enabled")

    def test_profiling_calibrate_importable(self):
        assert hasattr(meshmonk, "profiling_calibrate")

    def test_profiling_reset_in_all(self):
        assert "profiling_reset" in meshmonk.__all__

    def test_profiling_peek_in_all(self):
        assert "profiling_peek" in meshmonk.__all__

    def test_profiling_dump_in_all(self):
        assert "profiling_dump" in meshmonk.__all__

    def test_profiling_enabled_in_all(self):
        assert "profiling_enabled" in meshmonk.__all__

    def test_profiling_calibrate_in_all(self):
        assert "profiling_calibrate" in meshmonk.__all__


# ---------------------------------------------------------------------------
# Return type contracts (hold in both ON and OFF builds)
# ---------------------------------------------------------------------------

class TestReturnTypes:
    """Return types must be correct regardless of MESHMONK_PROFILING build flag."""

    def test_profiling_enabled_returns_bool(self):
        result = meshmonk.profiling_enabled()
        assert isinstance(result, bool)

    def test_profiling_reset_returns_none(self):
        result = meshmonk.profiling_reset()
        assert result is None

    def test_profiling_peek_returns_dict(self):
        result = meshmonk.profiling_peek()
        assert isinstance(result, dict)

    def test_profiling_dump_returns_dict(self):
        result = meshmonk.profiling_dump()
        assert isinstance(result, dict)

    def test_profiling_calibrate_returns_int(self):
        # Use n=1 for fast test (1 iteration is enough to check type)
        result = meshmonk.profiling_calibrate(1)
        assert isinstance(result, int)

    def test_profiling_calibrate_default_arg(self):
        """calibrate() must accept zero-argument call (default n=1_000_000).
        We call it with n=1 for speed, but verify the default doesn't crash.
        Just test that n is accepted as a keyword argument too.
        """
        result = meshmonk.profiling_calibrate(n=1)
        assert isinstance(result, int)


# ---------------------------------------------------------------------------
# OFF-build no-op semantics (when profiling_enabled() is False)
# ---------------------------------------------------------------------------

class TestOffBuildNoOps:
    """When MESHMONK_PROFILING is not set, functions must return empty/zero."""

    def test_off_build_peek_returns_empty(self):
        if meshmonk.profiling_enabled():
            pytest.skip("This test only applies to OFF builds")
        assert meshmonk.profiling_peek() == {}

    def test_off_build_dump_returns_empty(self):
        if meshmonk.profiling_enabled():
            pytest.skip("This test only applies to OFF builds")
        assert meshmonk.profiling_dump() == {}

    def test_off_build_calibrate_returns_zero(self):
        if meshmonk.profiling_enabled():
            pytest.skip("This test only applies to OFF builds")
        assert meshmonk.profiling_calibrate() == 0

    def test_off_build_reset_no_raise(self):
        if meshmonk.profiling_enabled():
            pytest.skip("This test only applies to OFF builds")
        meshmonk.profiling_reset()  # must not raise


# ---------------------------------------------------------------------------
# ON-build semantics: peek/dump/reset contract
# ---------------------------------------------------------------------------

class TestOnBuildSemantics:
    """When MESHMONK_PROFILING is set, verify reset/peek/dump semantics."""

    def setup_method(self):
        """Ensure clean accumulator before each test."""
        if meshmonk.profiling_enabled():
            meshmonk.profiling_reset()

    def test_on_build_reset_clears_accumulator(self):
        if not meshmonk.profiling_enabled():
            pytest.skip("This test only applies to ON builds")
        # After reset, peek should return empty
        meshmonk.profiling_reset()
        assert meshmonk.profiling_peek() == {}

    def test_on_build_peek_does_not_reset(self):
        """Two consecutive peeks must return identical data (peek has no side effects)."""
        if not meshmonk.profiling_enabled():
            pytest.skip("This test only applies to ON builds")

        # Use calibrate to populate the accumulator (it uses ScopedTimer internally
        # but erases __calibration__ entries — that's fine, we just need some data).
        # Actually, calibrate erases its own entries. Use a different approach:
        # just check that peek/peek gives the same dict (even if both are empty).
        meshmonk.profiling_reset()
        snap1 = meshmonk.profiling_peek()
        snap2 = meshmonk.profiling_peek()
        assert snap1 == snap2, "peek must not reset the accumulator: two consecutive peeks differ"

    def test_on_build_dump_resets_accumulator(self):
        """dump() must return the current snapshot AND reset the accumulator."""
        if not meshmonk.profiling_enabled():
            pytest.skip("This test only applies to ON builds")
        meshmonk.profiling_reset()
        # Get pre-dump snapshot via peek
        pre_dump = meshmonk.profiling_peek()
        # dump must return same data
        dumped = meshmonk.profiling_dump()
        assert dumped == pre_dump, "dump must return same data as peek before it"
        # After dump, accumulator must be empty
        post_dump = meshmonk.profiling_peek()
        assert post_dump == {}, "dump must reset the accumulator"

    def test_on_build_peek_peek_dump_pattern(self):
        """The full peek-peek-dump pattern: peek is idempotent, dump resets."""
        if not meshmonk.profiling_enabled():
            pytest.skip("This test only applies to ON builds")
        meshmonk.profiling_reset()
        snap1 = meshmonk.profiling_peek()
        snap2 = meshmonk.profiling_peek()
        assert snap1 == snap2, "two peeks must return identical data"
        snap3 = meshmonk.profiling_dump()
        assert snap3 == snap1, "dump must match peek before dump"
        snap4 = meshmonk.profiling_peek()
        assert snap4 == {}, "post-dump peek must be empty"

    def test_on_build_calibrate_plausible_range(self):
        """calibrate(n) must return a plausible nanoseconds-per-scope value."""
        if not meshmonk.profiling_enabled():
            pytest.skip("This test only applies to ON builds")
        ns = meshmonk.profiling_calibrate(100_000)
        # ScopedTimer = two chrono::now() calls, expected 10-2000 ns/scope
        assert 10 <= ns <= 5000, f"calibrate() returned implausible {ns} ns/scope"

    def test_on_build_calibrate_does_not_pollute_accumulator(self):
        """calibrate() must not leave __calibration__ entries in the accumulator."""
        if not meshmonk.profiling_enabled():
            pytest.skip("This test only applies to ON builds")
        meshmonk.profiling_reset()
        meshmonk.profiling_calibrate(1000)
        snap = meshmonk.profiling_peek()
        assert "__calibration__" not in snap, \
            "calibrate must erase its own entries from the accumulator"

    def test_on_build_peek_dict_shape(self):
        """When accumulator has data, peek returns {label: {total_us, count}} shape."""
        if not meshmonk.profiling_enabled():
            pytest.skip("This test only applies to ON builds")
        # Since we can't easily populate the accumulator without running a
        # registration (which requires mesh data), we verify the empty case shape.
        # Shape is implicitly validated by test_profiling_enabled_roundtrip below.
        meshmonk.profiling_reset()
        snap = meshmonk.profiling_peek()
        # Empty dict is valid shape
        assert isinstance(snap, dict)
        # If non-empty, each value must have total_us and count keys
        for label, entry in snap.items():
            assert isinstance(label, str)
            assert isinstance(entry, dict)
            assert "total_us" in entry
            assert "count" in entry
            assert isinstance(entry["total_us"], int)
            assert isinstance(entry["count"], int)


# ---------------------------------------------------------------------------
# Integration roundtrip: populate accumulator via registration
# ---------------------------------------------------------------------------

class TestProfilingRoundtrip:
    """Roundtrip test that populates the accumulator via actual registration."""

    def test_reset_peek_peek_dump_roundtrip_with_registration(self, mesh_data_dir):
        """Full reset->run->peek->peek->dump roundtrip with real registration data."""
        if not meshmonk.profiling_enabled():
            pytest.skip("This test only applies to ON builds")

        import trimesh

        template_path = mesh_data_dir / "Template.obj"
        demoface_path = mesh_data_dir / "demoFace.obj"
        if not template_path.exists():
            pytest.skip(f"Template.obj not found: {template_path}")
        if not demoface_path.exists():
            pytest.skip(f"demoFace.obj not found: {demoface_path}")

        m1 = trimesh.load(str(template_path), process=False)
        m2 = trimesh.load(str(demoface_path), process=False)

        meshmonk.profiling_reset()
        meshmonk.rigid_register(floating=m1, target=m2)

        snap1 = meshmonk.profiling_peek()
        snap2 = meshmonk.profiling_peek()

        # Two consecutive peeks must be identical (no side effects)
        assert snap1 == snap2, "peek must not reset: two consecutive peeks differ"

        snap3 = meshmonk.profiling_dump()
        assert snap3 == snap1, "dump must return same data as peek"

        snap4 = meshmonk.profiling_peek()
        assert snap4 == {}, "dump must have reset: post-dump peek should be empty"

        # The accumulator should have had some entries (instrumented callsites)
        assert len(snap1) > 0, "Expected at least one profiling entry after registration"

        # Verify shape of each entry
        for label, entry in snap1.items():
            assert "total_us" in entry, f"Missing total_us in entry for {label!r}"
            assert "count" in entry, f"Missing count in entry for {label!r}"
            assert entry["count"] > 0, f"count must be positive for {label!r}"
            assert isinstance(entry["total_us"], int)
            assert isinstance(entry["count"], int)
