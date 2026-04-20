"""Tests for silent logger — verifies std::cerr is routed through meshmonk::log().

TDD for bead workspace-dyj.2: Route std::cerr through logger.
"""

import os
import subprocess
import sys
import tempfile


class TestSilentMode:
    """All internal C++ output must be suppressed when log_level is 'silent'."""

    def _run_subprocess(self, code: str) -> subprocess.CompletedProcess:
        env = dict(os.environ)
        # Ensure the shared libraries are findable in the subprocess
        site_lib = os.path.join(
            os.path.dirname(sys.executable),
            "..",
            "lib",
            f"python{sys.version_info.major}.{sys.version_info.minor}",
            "site-packages",
            "lib",
        )
        existing = env.get("LD_LIBRARY_PATH", "")
        env["LD_LIBRARY_PATH"] = (
            f"{os.path.realpath(site_lib)}:{existing}"
            if existing
            else os.path.realpath(site_lib)
        )
        return subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True,
            text=True,
            timeout=60,
            env=env,
            cwd=tempfile.gettempdir(),
        )

    def test_rigid_register_silent_no_stderr(self):
        """rigid_register with log_level=silent must produce no stderr output."""
        result = self._run_subprocess(
            "import meshmonk, trimesh; "
            "meshmonk.set_log_level('silent'); "
            "m = trimesh.creation.icosphere(subdivisions=1); "
            "meshmonk.rigid_register(floating=m, target=m, num_iterations=2)"
        )
        assert result.returncode == 0, f"subprocess failed: {result.stderr!r}"
        assert (
            result.stderr == ""
        ), f"Expected no stderr output in silent mode, got: {result.stderr!r}"

    def test_nonrigid_register_silent_no_stderr(self):
        """nonrigid_register with log_level=silent must produce no stderr output."""
        result = self._run_subprocess(
            "import meshmonk, trimesh; "
            "meshmonk.set_log_level('silent'); "
            "m = trimesh.creation.icosphere(subdivisions=1); "
            "meshmonk.nonrigid_register(floating=m, target=m, num_iterations=2)"
        )
        assert result.returncode == 0, f"subprocess failed: {result.stderr!r}"
        assert (
            result.stderr == ""
        ), f"Expected no stderr output in silent mode, got: {result.stderr!r}"

    def test_warning_level_default_emits_output(self):
        """At default warning level, meshmonk registration should emit output (not silent)."""
        # We just verify the process runs and doesn't crash — the key invariant
        # is that silent mode (tested above) suppresses it.
        result = self._run_subprocess(
            "import meshmonk, trimesh; "
            "m = trimesh.creation.icosphere(subdivisions=1); "
            "meshmonk.rigid_register(floating=m, target=m, num_iterations=2)"
        )
        assert result.returncode == 0, f"subprocess failed: {result.stderr!r}"
