"""Tests for silent logger — verifies std::cerr is routed through meshmonk::log().

TDD for bead workspace-dyj.2: Route std::cerr through logger.
"""

import os
import subprocess
import sys
import tempfile


# compute_inlier_weights with flipped normals logs a [WARNING] at the default level
# (InlierDetector.cpp); registration on identical icospheres logs nothing.
_FLIPPED = (
    "import numpy as np; "
    "V = np.asarray(m.vertices, np.float32); "
    "N = np.asarray(m.vertex_normals, np.float32); "
    "meshmonk.compute_inlier_weights("
    "np.hstack([V, N]), np.hstack([V, -N]), np.ones(len(V), np.float32))"
)


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

    def test_compute_inlier_weights_silent_after_rigid_register(self):
        """Silent mode suppresses compute_inlier_weights warnings after rigid registration."""
        result = self._run_subprocess(
            "import meshmonk, trimesh; "
            "meshmonk.set_log_level('silent'); "
            "m = trimesh.creation.icosphere(subdivisions=1); "
            "meshmonk.rigid_register(floating=m, target=m, num_iterations=2); "
            + _FLIPPED
        )
        assert result.returncode == 0, f"subprocess failed: {result.stderr!r}"
        assert (
            result.stderr == ""
        ), f"Expected no stderr output in silent mode, got: {result.stderr!r}"

    def test_compute_inlier_weights_silent_after_nonrigid_register(self):
        """Silent mode suppresses compute_inlier_weights warnings after nonrigid registration."""
        result = self._run_subprocess(
            "import meshmonk, trimesh; "
            "meshmonk.set_log_level('silent'); "
            "m = trimesh.creation.icosphere(subdivisions=1); "
            "meshmonk.nonrigid_register(floating=m, target=m, num_iterations=2); "
            + _FLIPPED
        )
        assert result.returncode == 0, f"subprocess failed: {result.stderr!r}"
        assert (
            result.stderr == ""
        ), f"Expected no stderr output in silent mode, got: {result.stderr!r}"

    def test_compute_inlier_weights_default_warning_level_emits_warning(self):
        """At the default warning level, compute_inlier_weights emits a warning."""
        result = self._run_subprocess(
            "import meshmonk, trimesh; "
            "m = trimesh.creation.icosphere(subdivisions=1); "
            "meshmonk.rigid_register(floating=m, target=m, num_iterations=2); "
            + _FLIPPED
        )
        assert result.returncode == 0, f"subprocess failed: {result.stderr!r}"
        assert "[WARNING]" in result.stderr, f"expected a warning: {result.stderr!r}"
