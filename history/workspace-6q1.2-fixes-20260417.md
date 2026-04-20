# Fixes: workspace-6q1.2
Date: 2026-04-17
Review passes: 1

## Critical
(none)

## Important
(none)

## Minor

- **library/src/meshmonk.cpp:54-73** — `run_final_inlier_pass` uses raw `new`/`delete` for the correspondence filter. The same cleanup (raw → `std::unique_ptr`) was applied to `RigidRegistration.cpp` and `NonrigidRegistration.cpp` as part of this change but was not applied to this new helper function. Not a bug (delete is always reached before return), but inconsistent with the cleanup intent. Fix: `auto corr_filter = std::unique_ptr<registration::BaseCorrespondenceFilter>{...}` and remove the explicit `delete corr_filter`.

- **library/src/meshmonk.cpp:111-112 / 182-183 / 263-264** — `FacesMat = Eigen::Matrix<int, Dynamic, 3>` has a compile-time column count of 3, so `.cols() != 3` is always false. The DegenerateInput check for column counts is dead code. Low impact (it builds and never fires incorrectly), but it gives a false sense of validation. In v0.2, consider a static_assert or note that this guard is intentionally defensive against future API changes. Not worth fixing now but worth a `// cols always == 3 (fixed at compile time); kept for documentation` comment.

- **library/src/meshmonk.cpp:157** — `result.iterations_run = params.num_iterations` is set to the *requested* number, not the *actual* iterations run (which is always `params.num_iterations` for fixed-iteration ICP). This is semantically correct for v0.1 since there is no early exit, but the field name implies actual iterations. Add a comment: `// v0.1: fixed-iteration ICP; actual == requested`. Same pattern applies to NonrigidResult (line 242).

- **library/include/meshmonk/meshmonk.hpp** — `compute_correspondences` and related primitives are not `[[nodiscard]]`. The high-level functions are correctly marked; omitting it from value-returning primitives is a missed opportunity for the compiler to catch ignored results. Consider adding `[[nodiscard]]` to all free functions that return meaningful values.

## ADR Updates
No ADR changes needed. All implementation choices are consistent with ADR-001 D6 (params structs, tl::expected, value returns, internal classes untouched, boundary error detection).

## Discarded

- **InsufficientInliers absent in rigid_registration**: The bead spec says this should be checked "after InlierDetector runs (in the first iteration)". However, ADR D6 explicitly prohibits rewriting internal classes. Raising `InsufficientInliers` for rigid would require either a post-registration inlier pass (extra computation on the already-converged mesh — not specified for rigid) or introspecting inside `RigidRegistration::update()` (requires rewriting load-bearing code). The nonrigid and pyramid paths correctly do a post-registration pass and check there. For rigid, the check is structurally absent — this is a known gap, not a code defect. Filed as discarded because the fix requires a design-level decision about whether rigid should also do a post-registration inlier pass.

- **RigidTransform::inverse() math concern**: Reviewed `Rt / s2` in transform.hpp. For pure rotation s=1 → `R^T/1 = R^T` ✓. For similarity with `R = sR_true`, `Rt = s*R_true^T`, `s2 = s^2`, so `Rt/s2 = (s*R_true^T)/s^2 = R_true^T/s = (1/s)R_true^T` ✓. The math is correct; discarded.

- **FacesMat cols check produces dead code**: Listed under Minor above; not discarded — just low priority.

- **compute_normals type compatibility**: `meshmonk::Vec3Mat` and `registration::Vec3Mat` both resolve to `Eigen::Matrix<float, Dynamic, 3>`. No type mismatch at the call site. Discarded.
