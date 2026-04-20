# Fixes: workspace-6q1.3
Date: 2026-04-17
Review passes: 1

## Critical
(none)

## Important
(none)

## Minor
- **meshmonk/__init__.py:264-270** — `_check_normals_and_warn` is dead code. It defines a function that does nothing (returns `features` unchanged) and is never called. The actual all-zero normals check and recompute is correctly inlined in each Pattern B branch. Safe to delete the stub to avoid confusion.
- **meshmonk/__init__.py:432-433** — Mixed-input guard. The condition `if floating is not None or target is not None` means if only one of `floating`/`target` is passed (e.g., `rigid_register(floating=mesh)`) the code enters Pattern A and then calls `_mesh_to_arrays(None, ...)`, which blows up with an AttributeError at `None.vertices`. The error message will be cryptic. A cheap guard at the top of each registration function would give a clear message: `if (floating is None) != (target is None): raise ValueError("'floating' and 'target' must both be provided or both omitted")`. Applies to all three register functions.
- **meshmonk/__init__.py (Pattern B)** — Target normals not checked. Pattern B checks `feat_float[:, 3:]` for all-zero and warns/recomputes, but does not check `feat_target[:, 3:]`. A user passing all-zero target normals proceeds silently with degenerate 6D correspondence search. Consistent with spec (spec only specifies floating check) but the symmetry gap is worth noting. Low priority — can be deferred to v0.2.

## ADR Updates
No ADR changes needed. All implementation decisions match ADR-001 (D2 nanobind, D6 API design) and the v0.1 design doc.

## Discarded
- **`set_log_level` missing GIL call_guard**: The bead spec says GIL release is "MANDATORY for all registration functions AND all primitives." `set_log_level` does a single atomic write — adding `nb::call_guard<nb::gil_scoped_release>()` would be technically correct but offers zero practical benefit. Discarded as noise; not worth the churn.
- **`inverse()` scale formula in transform.hpp**: `s2 = R.col(0).squaredNorm()` computes s² correctly for the similarity case. However `inv.topLeftCorner<3,3>() = Rt / s2` yields `R^T / s²` where for a similarity transform `sR` the correct inverse rotation block is `(1/s) R^T = R^T / s`. This looks like a bug (divides by s² instead of s), but it is in the C++ library itself (bead 2 scope), not the bindings. Out of scope for bead 3; filed separately if confirmed.
- **`NonrigidParams` missing `use_scaling` field in kwarg flattening**: `use_scaling` is only on `RigidParams`, not `NonrigidParams` or `PyramidParams`, which is correct per the design and params.hpp. Not a missing binding.
