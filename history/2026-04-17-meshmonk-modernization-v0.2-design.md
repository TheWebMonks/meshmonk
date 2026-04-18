# Design: MeshMonk v0.2 — Team-ready

**Date:** 2026-04-17
**Status:** Draft
**Parent:** [Overview](./2026-04-17-meshmonk-modernization-design.md)
**Decisions:** [ADR-001](../docs/decisions/ADR-001-meshmonk-modernization.md)

---

## Scope

**v0.2 — Team-ready:** ship type stubs. README rewrite. Migration guide from MATLAB (side-by-side snippets) at `docs/migration-from-matlab.md`. Windows CI if feasible. cibuildwheel dry-run. Add diagnostic fields (`converged`, `fitness`, `inlier_rmse`, `num_inliers`) to result structs once convergence criteria are pinned down. PR `cli` → `master`; rename `master` → `main`.

## Deferred from v0.1 review

Items identified during the v0.1 code review (2026-04-18) that were out of scope for v0.1 but should land in v0.2. See `history/meshmonk-modernization-v0.1-fixes-20260417.md` and `history/v01-review-fixes-20260418.md` for full context.

### API cleanup

- **Extract `_prepare_arrays()` helper** — `rigid_register`, `nonrigid_register`, and `pyramid_register` each contain ~75 lines of identical Pattern A/B validation and array-preparation logic (~225 lines total). Extract into a shared `_prepare_arrays(floating, target, *, floating_features, ..., normals, compute_normals_flag) -> (feat_float, feat_target, faces_float, faces_target, flags_float, flags_target)` helper. Each public function then calls `_prepare_arrays()` followed by its specific `_apply_*_kwargs` and `_*_registration`.

- **`rigid_params` kwarg for nonrigid/pyramid** — Allow explicit control of the rigid pre-alignment stage in `nonrigid_register()` and `pyramid_register()`. Requires C++ API changes: `nonrigid_registration()` and `pyramid_registration()` currently have no separate rigid params input.

### Error handling

- **Typed `MeshMonkError` exception class** — Currently `MeshMonkError` is bound to `std::runtime_error` in nanobind, catching ALL `std::runtime_error` from the entire call stack (OpenMesh, Eigen, internal classes). Callers must string-parse the error message to distinguish `DegenerateInput` from `InsufficientInliers`. Fix: create a dedicated C++ `MeshMonkError` class carrying a `RegistrationError` enum field, exposed as `error_code` in Python.

- **Route `std::cerr` through logger sink** — Internal legacy classes write warnings/errors directly to `std::cerr`, bypassing `set_log_level("silent")`. Users calling `pyramid_register()` on small meshes get Downsampler warnings to stderr that cannot be suppressed. Route through the existing logger sink so `set_log_level` controls all output.

### CLI

- **`demo --download` URLs** — The `demo --download` subcommand is a stub that prints "TODO: set download URLs" and exits with code 0 (success). Either implement actual download URLs or exit with code 1 and a clear message.
