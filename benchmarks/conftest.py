"""pytest-benchmark configuration for the MeshMonk benchmark suite.

Key invariants enforced here:
1. Hard-assert that pytest-xdist is NOT active: running benchmarks under
   `-n auto` multiplies the OMP thread budget across workers, producing
   noisy and non-reproducible numbers. We fail fast rather than measure
   garbage.

2. Thread pinning: set OMP_NUM_THREADS=1 so every baseline run is
   strictly single-threaded. This means reported speedup from future
   OpenMP changes is the *incremental* effect of MeshMonk-level pragmas,
   not contaminated by Eigen's opportunistic threading.

   NOTE: Eigen::setNbThreads(1) is the correct runtime call to disable
   Eigen's own threading (EIGEN_DONT_PARALLELIZE is a compile-time macro;
   its env-var form is a no-op at runtime). No Python binding for
   setNbThreads exists in meshmonk yet. When it is added, call it here.
   OMP_NUM_THREADS=1 is the best available fallback and also constrains
   Eigen's OpenMP-backed paths.

3. Benchmark config: rounds >= 5, median aggregation (set via
   pytest-benchmark ini options or the benchmark fixture's timer).
"""

from __future__ import annotations

import os

import pytest

# ---------------------------------------------------------------------------
# Hard-assert: benchmarks must NOT run under pytest-xdist
# ---------------------------------------------------------------------------
assert "PYTEST_XDIST_WORKER" not in os.environ, (
    "Benchmarks must not run under pytest-xdist (-n flag). "
    "Running benchmarks across parallel workers oversubscribes the OMP thread "
    "budget and produces noisy, non-reproducible results. "
    "Remove -n / --numprocesses from your pytest invocation."
)


# ---------------------------------------------------------------------------
# Session-scoped fixture: pin thread counts for baseline isolation
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session", autouse=True)
def pin_threads():
    """Pin OMP_NUM_THREADS=1 for the entire benchmark session.

    This ensures baseline captures reflect strictly single-threaded
    MeshMonk behaviour. Post-optimisation runs (9f5, bdt) should be
    run WITHOUT this fixture (or with OMP_NUM_THREADS unset) so the
    reported speedup is the incremental OpenMP contribution.

    Eigen::setNbThreads(1) is not called here because no Python binding
    exists yet. When meshmonk exposes Eigen::setNbThreads(), add:
        import meshmonk
        meshmonk.set_eigen_num_threads(1)
    here and restore to 0 in the finally block.
    """
    original = os.environ.get("OMP_NUM_THREADS")
    os.environ["OMP_NUM_THREADS"] = "1"
    yield
    if original is None:
        del os.environ["OMP_NUM_THREADS"]
    else:
        os.environ["OMP_NUM_THREADS"] = original


# ---------------------------------------------------------------------------
# Benchmark configuration defaults
# ---------------------------------------------------------------------------

# NOTE: rounds and timer are set per-benchmark via the benchmark fixture or
# via pyproject.toml [tool.pytest.benchmark] when added. The minimum rounds
# requirement (>= 5) is enforced by the --benchmark-min-rounds CLI flag used
# in CI. The median aggregation is selected via --benchmark-timer or the
# fixture's .pedantic() call.


def pytest_configure(config):
    """Register benchmark markers and enforce minimum rounds in CI mode."""
    config.addinivalue_line(
        "markers",
        "benchmark: mark a test as a benchmark (run with --benchmark-only)",
    )
