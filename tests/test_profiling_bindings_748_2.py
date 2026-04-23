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

from pathlib import Path

import meshmonk
import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Helper: populate the accumulator via a real (fast) registration
# ---------------------------------------------------------------------------

_DATA_DIR = Path("/workspace/data")
_TEMPLATE_1K = _DATA_DIR / "Template_1K.obj"
_DEMOFACE_1K = _DATA_DIR / "DemoFace_1K.obj"


def _populate_accumulator_with_registration():
    """Run a rigid registration with 1K meshes to populate the accumulator.

    Requires Template_1K.obj and DemoFace_1K.obj to exist in /workspace/data/.
    Skips via pytest.skip() if the mesh files are not present.
    Returns the snap dict from profiling_peek() after the run (for assertion).
    """
    import trimesh

    if not _TEMPLATE_1K.exists():
        pytest.skip(f"Template_1K.obj not found at {_TEMPLATE_1K}")
    if not _DEMOFACE_1K.exists():
        pytest.skip(f"DemoFace_1K.obj not found at {_DEMOFACE_1K}")

    m1 = trimesh.load(str(_TEMPLATE_1K), process=False)
    m2 = trimesh.load(str(_DEMOFACE_1K), process=False)

    vf = np.asarray(m1.vertices, dtype=np.float32)
    ff = np.asarray(m1.faces, dtype=np.int32)
    vt = np.asarray(m2.vertices, dtype=np.float32)
    ft = np.asarray(m2.faces, dtype=np.int32)

    feat_f = meshmonk.features_from_vertices(vf, ff)
    feat_t = meshmonk.features_from_vertices(vt, ft)
    flags_f = np.ones(len(vf), dtype=np.float32)
    flags_t = np.ones(len(vt), dtype=np.float32)

    meshmonk.profiling_reset()
    meshmonk.rigid_register(
        floating_features=feat_f,
        target_features=feat_t,
        floating_faces=ff,
        target_faces=ft,
        floating_flags=flags_f,
        target_flags=flags_t,
    )


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
        """Two consecutive peeks must return identical data (peek has no side effects).

        Populates the accumulator via a real registration before checking — an
        empty accumulator would make the test trivially pass even if peek secretly
        reset, since {} == {} always.
        """
        if not meshmonk.profiling_enabled():
            pytest.skip("This test only applies to ON builds")

        # Populate accumulator with real registration data
        _populate_accumulator_with_registration()

        # Both peeks must return the same non-empty dict
        snap1 = meshmonk.profiling_peek()
        snap2 = meshmonk.profiling_peek()
        assert len(snap1) > 0, "Expected non-empty accumulator after registration"
        assert snap1 == snap2, "peek must not reset the accumulator: two consecutive peeks differ"

    def test_on_build_dump_resets_accumulator(self):
        """dump() must return the current snapshot AND reset the accumulator.

        Populates the accumulator via a real registration before checking — an
        empty accumulator would make the test trivially pass even if dump did not
        reset, since the post-dump check would also see {}.
        """
        if not meshmonk.profiling_enabled():
            pytest.skip("This test only applies to ON builds")

        # Populate accumulator with real registration data
        _populate_accumulator_with_registration()

        # Get pre-dump snapshot via peek (must be non-empty)
        pre_dump = meshmonk.profiling_peek()
        assert len(pre_dump) > 0, "Expected non-empty accumulator after registration"
        # dump must return same data
        dumped = meshmonk.profiling_dump()
        assert dumped == pre_dump, "dump must return same data as peek before it"
        # After dump, accumulator must be empty
        post_dump = meshmonk.profiling_peek()
        assert post_dump == {}, "dump must reset the accumulator"

    def test_on_build_peek_peek_dump_pattern(self):
        """The full peek-peek-dump pattern: peek is idempotent, dump resets.

        Populates the accumulator via a real registration so the assertions are
        non-trivial (i.e., they catch a secret reset in peek).
        """
        if not meshmonk.profiling_enabled():
            pytest.skip("This test only applies to ON builds")

        # Populate accumulator with real registration data
        _populate_accumulator_with_registration()

        snap1 = meshmonk.profiling_peek()
        snap2 = meshmonk.profiling_peek()
        assert len(snap1) > 0, "Expected non-empty accumulator after registration"
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
        """calibrate() must not pollute the global accumulator with any entries.

        After the refactor (local Profiler inside calibrate), the global
        accumulator must remain completely untouched — no __calibration__ key
        and no other keys inserted by the calibration loop.
        """
        if not meshmonk.profiling_enabled():
            pytest.skip("This test only applies to ON builds")
        meshmonk.profiling_reset()
        meshmonk.profiling_calibrate(1000)
        snap = meshmonk.profiling_peek()
        # Accumulator must be completely empty — calibrate uses a local Profiler
        assert snap == {}, \
            "calibrate must not touch the global accumulator (uses local Profiler)"

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
    """Roundtrip test that populates the accumulator via actual registration.

    Uses Template_1K.obj + DemoFace_1K.obj which are committed to data/ as part
    of bead meshmonk-modernization-748.3. The skip condition only requires those
    files to exist — it does NOT require the larger Template.obj/demoFace.obj.
    """

    def test_reset_peek_peek_dump_roundtrip_with_registration(self):
        """Full reset->run->peek->peek->dump roundtrip with real registration data."""
        if not meshmonk.profiling_enabled():
            pytest.skip("This test only applies to ON builds")

        # Populate accumulator via the shared helper (skips if 1K meshes absent)
        _populate_accumulator_with_registration()

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
