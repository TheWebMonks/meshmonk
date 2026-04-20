# Fixes: meshmonk-modernization-v0.1
Date: 2026-04-17
Review passes: 2 (architecture + implementation per pass)

## Critical

- **library/include/meshmonk/transform.hpp:42-46** — `RigidTransform::inverse()` divides by `s^2` instead of `s`. Correct inverse of `[sR|t]` is `[(1/s)R^T | -(1/s)R^T*t]`. Current code: `Rt / s2` produces `R^T/s^2`. Fix: `float s = std::sqrt(s2); inv.topLeftCorner<3,3>() = Rt / s; inv.topRightCorner<3,1>() = -(Rt * t) / s;`. Only affects `use_scaling=true` (s=1 case is unaffected).

- **library/CMakeLists.txt:16** — `meshmonk_lib` is `SHARED` but no `install(TARGETS meshmonk_lib ...)` exists. pip-installed wheels contain `_meshmonk_core.so` that links against `libmeshmonk_lib.so` which is NOT in the wheel → runtime `ImportError`. Fix: change `add_library(meshmonk_lib SHARED ...)` to `add_library(meshmonk_lib STATIC ...)` so it's statically linked into `_meshmonk_core.so`.

- **library/src/meshmonk.cpp** — Division by zero when `num_iterations=1` in NonrigidRegistration. `NonrigidRegistration.cpp:53` computes `exp(log(end/start) / (numIterations - 1))`. Fix: add guard in `nonrigid_registration()` and `pyramid_registration()`: `if (params.num_iterations < 2) return tl::unexpected{RegistrationError::DegenerateInput};` or clamp to minimum 2.

## Important

- **library/src/meshmonk.cpp:54-73,369-383** — `run_final_inlier_pass()` and `compute_correspondences()` use raw `new`/`delete` for `BaseCorrespondenceFilter*`. Exception-unsafe. Fix: `auto filter = std::unique_ptr<registration::BaseCorrespondenceFilter>{...};` and remove `delete` calls.

- **library/include/meshmonk/transform.hpp:22-33** — `RigidTransform::apply()` uses `matrix.topLeftCorner<3,3>()` for normals, which is `sR` when `use_scaling=true`. Normals get scaled by `s` instead of just rotated. Fix: extract `float s = R.col(0).norm(); auto R_pure = R / s;` and apply `R_pure` to normals columns 3-5.

- **library/include/meshmonk/params.hpp:52** — `PyramidParams` inherits `flag_threshold=0.9` from `CorrespondenceParams`, but MATLAB pyramid uses `0.999`. Python wrapper patches this but C++ callers get wrong default. Fix: add `PyramidParams() { correspondences.flag_threshold = 0.999f; }` constructor.

- **meshmonk/__init__.py:432,499,567** — Pattern A crashes with `AttributeError` when only one of `floating`/`target` is provided. Fix: add guard `if (floating is None) != (target is None): raise ValueError("both floating and target must be provided")` before Pattern A entry.

- **meshmonk/__init__.py:273-371** — No validation of unrecognized kwargs. `inlier_kapaa=12.0` (typo) silently uses default. Fix: add known-kwargs set per function, raise `TypeError` for unrecognized keys.

- **library/src/meshmonk.cpp:175-186,258-270** — Empty faces (`rows()==0`) not guarded for nonrigid/pyramid. Fix: add `if (floating_faces.rows() == 0 || target_faces.rows() == 0) return tl::unexpected{RegistrationError::DegenerateInput};`.

- **library/src/meshmonk.cpp:130-158** — `InsufficientInliers` not checked after rigid registration (inconsistent with nonrigid/pyramid). Fix: add `run_final_inlier_pass()` and the `<4` check in `rigid_registration()`.

- **meshmonk/__init__.py:514-524** — Target features all-zero normals not auto-recomputed in Pattern B. Only floating is checked. Fix: apply the same all-zero check + auto-recompute to `feat_target` using `target_faces`.

- **library/src/meshmonk.cpp:113** — Missing `target_flags` all-zero guard. Fix: add `if (target_flags.maxCoeff() == 0.0f) return tl::unexpected{RegistrationError::DegenerateInput};`.

- **library/src/meshmonk.cpp:116-117** — Bounding box degeneracy check only on X-axis. Fix: check all 3 position columns (0, 1, 2) — degenerate if ALL three have zero range.

- **meshmonk/cli.py:162** — pyramid `--iterations` help text says "per layer" but value is total iterations. Fix: change to `help="Total iterations across all pyramid layers."`.

- **tests/test_cli.py:178** — `test_demo_runs_rigid_mode` is a false positive: `"--download" in output` always matches help text. Fix: remove the `"--download" in output` clause from the assertion.

## Minor

- **meshmonk/__init__.py:264-270** — Dead function `_check_normals_and_warn` (defined, never called, does nothing). Fix: delete it.

- **meshmonk/__init__.py** — `MeshMonkError` not re-exported from `meshmonk` package. Fix: add to imports from `_meshmonk_core` and include in module scope.

- **library/src/meshmonk.cpp:317-320** — `per_layer_iterations` approximation can diverge from `num_iterations` (e.g., 100/3 → 33*3=99). Fix: add comment documenting this is an approximation; consider distributing remainder to first layer.

## ADR Updates
- No ADR changes needed. All findings are implementation bugs or gaps, not design-level conflicts.

## Discarded

- **Triplicated validation logic** (Arch P1): Acceptable duplication. Design doc says "at TOP of each function". Extracting a helper is v0.2 refactor, not a fix.
- **scale_shift_mesh semantics mismatch** (Arch P1): Design doc Python sketch is illustrative. ScaleShifter works correctly for its actual (pyramid-internal) use.
- **Global include_directories()** (Arch P1): Works fine. Minor CMake hygiene, not a bug.
- **iterations_run == num_iterations** (Arch P1): Known, documented, becomes meaningful in v0.2 with convergence.
- **Missing rigid_params= kwarg** (Arch P1): v0.2 enhancement. Kwargs path is complete.
- **Missing __main__.py** (Arch P1): Nice-to-have, not a bug.
- **set_log_level GIL inconsistency** (Impl P1): Intentional — atomic store doesn't need GIL release.
- **exit(1) in internal code** (Arch P2): Real concern but in legacy internal code explicitly out of v0.1 modification scope. Document as known risk for v0.2.
- **std::cout debug spew** (Arch P2, Impl P2): In legacy internal classes. Modifying 5+ internal .cpp files is out of v0.1 scope. Document as v0.2 cleanup.
- **Wrong variable in PyramidNonrigidRegistration.cpp:54** (Impl P2): Real bug in legacy internal code. Out of v0.1 modification scope. File as separate issue.
- **Nightly golden tests only ubuntu** (Arch P2): Goldens are xfail stubs. Matters when goldens are captured (human task).
- **CI missing Ubuntu 24.04 + clang** (Arch P2): Minor CI gap, not blocking v0.1.
- **DecompositionFailed post-hoc only** (Arch P2): Known limitation. Internal code modification out of scope.
- **RigidTransform.matrix read-only** (Arch P2): v0.2 enhancement.
- **Pattern A ignores .flags** (Arch P2): v0.2 enhancement.
- **downsample_mesh ratio convention** (Arch P2): Documentation-only. Not a bug.
- **rigid_registration validates faces but doesn't use them** (Impl P2): API consistency design choice. Faces needed for future convergence check.
- **test_compute_rigid_transform_cube translation** (Impl P2): Test correctly scoped to what primitive can verify in single-shot.
- **nb::exception catches all runtime_errors** (Arch P2): MeshMonkError inherits RuntimeError. Users catch either. Minor semantic issue handled by re-export fix.
- **No Tier 1/2 nonrigid/pyramid tests** (Impl P1): Substantial new work. File as separate issue for v0.2 test expansion.
