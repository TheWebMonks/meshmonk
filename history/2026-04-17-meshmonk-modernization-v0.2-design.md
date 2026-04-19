# Design: MeshMonk v0.2 — Team-ready

**Date:** 2026-04-17
**Status:** Reviewed (3 passes, 2026-04-18)
**Parent:** [Overview](./2026-04-17-meshmonk-modernization-design.md)
**Decisions:** [ADR-001](../docs/decisions/ADR-001-meshmonk-modernization.md)

---

## Scope

**v0.2 — Team-ready:** ship type stubs. README rewrite. Migration guide from MATLAB (side-by-side snippets) at `docs/migration-from-matlab.md`. Windows CI assessment and implementation. cibuildwheel dry-run. PR `meshmonk-modernization` → `cli` → `master`; rename `master` → `main`. API cleanup and error handling improvements deferred from v0.1 review.

**Explicitly deferred to v0.3:**
- **Diagnostic fields** (`converged`, `fitness`, `inlier_rmse`, `num_inliers`) — prerequisites not met. No convergence criterion is defined in the legacy code (runs fixed iterations). `fitness` has no codebase definition. `num_inliers` requires a threshold on the continuous `InlierDetector` weights. These are research/specification tasks, not engineering tasks. Defer until concrete definitions are proposed and accepted.
- **`std::expected` migration** — gcc-11 and clang-15 in the current CI matrix do not support `<expected>`. Raising the compiler floor to gcc-13/clang-16 drops Ubuntu 22.04 coverage. `tl::expected` continues to work correctly. Defer until v0.3 when the compiler floor can be deliberately raised.
- **Pluggable logger sink** — `meshmonk::log()` writes to `std::cerr`, which is correct for v0.2 (the `std::cerr` routing work item routes _uncontrolled_ writes through the logger). A pluggable output backend (Python callback, file) is deferred to v0.3+.

## Implementation ordering

Items are grouped by priority. **Merge-preparation items** must land before the PR to `master` since they affect what visitors see.

### Priority 1 — Merge preparation

These items are prerequisites for the public-facing PR to master:

1. **README rewrite** (see [README scope](#readme-rewrite-scope) below)
2. **Migration guide** — `docs/migration-from-matlab.md` with side-by-side MATLAB → Python snippets
3. **Delete `cli/` directory** — ADR D4 (FIRM) says "Delete the entire `cli/` directory." `cli/CMakeLists.txt` and `cli.cpp` still exist. The root `CMakeLists.txt` no longer references them. Remove as dead code.
4. **Version bump** — bump `pyproject.toml` and `meshmonk/__init__.py` from `"0.1.0"` to `"0.2.0"` as the final commit before the PR to `master`
5. **`demo --download` fix** — change exit code from 0 (success) to 1 (error) with a clear message: "Demo meshes are bundled in the package at `meshmonk/data/`. Download URLs not yet configured." Actual download URLs deferred to v0.3 when a GitHub Release exists to host assets.
6. **Typed `MeshMonkError` exception class** (see [Error handling](#error-handling) below)
7. **Branch operations** (see [Branch operations checklist](#branch-operations-checklist) below)

### Priority 2 — API and code quality improvements

These improve the codebase but do not block the merge:

8. **Extract `_prepare_arrays()` helper** (see [API cleanup](#api-cleanup) below)
9. **`rigid_params` kwarg for nonrigid/pyramid** (see [API cleanup](#api-cleanup) below)
10. **Route `std::cerr` through logger sink** (see [Error handling](#error-handling) below)

### Priority 3 — Infrastructure and developer experience

Independent of the merge, can be done in parallel:

11. **Type stubs** (see [Type stubs](#type-stubs) below)
12. **`py.typed` marker** — add empty `meshmonk/py.typed` to signal PEP 561 compliance
13. **Windows CI** (see [Windows CI](#windows-ci) below)
14. **cibuildwheel dry-run** (see [cibuildwheel dry-run](#cibuildwheel-dry-run) below)

## Deferred from v0.1 review

Items identified during the v0.1 code review (2026-04-18) that were out of scope for v0.1 but should land in v0.2. See `history/meshmonk-modernization-v0.1-fixes-20260417.md` and `history/v01-review-fixes-20260418.md` for full context.

### API cleanup

- **Extract `_prepare_arrays()` helper** — `rigid_register`, `nonrigid_register`, and `pyramid_register` each contain ~75 lines of identical Pattern A/B validation and array-preparation logic (~225 lines total). Extract into a shared helper. Each public function then calls `_prepare_arrays()` followed by its specific `_apply_*_kwargs` and C++ call.

  Refined signature (note: `normals` applies to floating mesh only, not both):
  ```python
  _prepare_arrays(
      floating, target, *,
      floating_features=None, target_features=None,
      floating_normals_override=None,  # applies to floating only
      compute_normals=False,
      floating_flags=None, target_flags=None,
  ) -> (feat_float, feat_target, faces_float, faces_target, flags_float, flags_target)
  ```

  `_prepare_arrays()` handles Python-level dispatch (Pattern A/B) and type coercion only. C++ `DegenerateInput` pre-checks remain in `meshmonk.cpp` to protect direct C++ callers.

- **`rigid_params` kwarg for nonrigid/pyramid** — Allow explicit control of the rigid pre-alignment stage in `nonrigid_register()` and `pyramid_register()`. **Implementation: Python-level composition** — run `rigid_registration` first, then feed its output into `nonrigid_registration`/`pyramid_registration`. This avoids touching internal C++ classes (`NonrigidRegistration` and `PyramidNonrigidRegistration` have no rigid pre-alignment stage — the MATLAB workflow calls rigid separately). No C++ API changes required.

### Error handling

- **Typed `MeshMonkError` exception class** — Currently `MeshMonkError` is bound to `std::runtime_error` in nanobind, catching ALL `std::runtime_error` from the entire call stack (OpenMesh, Eigen, internal classes). Callers must string-parse the error message to distinguish `DegenerateInput` from `InsufficientInliers`. Fix: create a dedicated C++ `MeshMonkError` class (`class MeshMonkError : public std::runtime_error` with `RegistrationError code` member), update `unwrap_expected` to throw this class, and bind with `nb::exception<MeshMonkError>(...)`.

  **Behavioral change (intentional):** After this fix, `MeshMonkError` will only catch throws of the custom `MeshMonkError` class. `std::runtime_error` from OpenMesh, Eigen, or other dependencies will surface as Python `RuntimeError`. This is better behavior — users can distinguish library errors from unexpected crashes — but is a breaking change for callers who catch `MeshMonkError` expecting all C++ errors. Acceptable at pre-1.0. The migration guide should mention this.

- **Route `std::cerr` through logger sink** — Internal legacy classes write warnings/errors directly to `std::cerr`, bypassing `set_log_level("silent")`. Route through the existing `meshmonk::log()` sink so `set_log_level` controls all output.

  **Scope:** 20 active `std::cerr` writes across 7 internal source files (`InlierDetector.cpp`, `Downsampler.cpp`, `ViscoElasticTransformer.cpp`, `BaseCorrespondenceFilter.cpp`, `helper_functions.cpp`, `ScaleShifter.cpp`, `RigidTransformer.cpp`). Implementation is a mechanical replacement (`std::cerr << msg` → `meshmonk::log(LogLevel::Warning, msg)`) that does not change control flow. Respects ADR D6 constraint: "internal algorithmic classes are NOT rewritten."

  **Note:** `meshmonk::log()` itself writes to `std::cerr` — this is intentional. The value is centralized control via `set_log_level`, not redirecting output away from stderr.

### CLI

- **`demo --download` fix** — The `demo --download` subcommand is a stub that prints "TODO: set download URLs" and exits with code 0 (success). Fix: exit with code 1 and a clear message. Actual download URLs deferred to v0.3 when demo meshes can be hosted as GitHub Release assets.

## v0.1-deferred items reconciliation

Three items from `history/meshmonk-modernization-v0.1-fixes-20260417.md` "Discarded" section were explicitly deferred to v0.2 but not captured in the original v0.2 design:

| Item | Disposition |
|------|-------------|
| **"Wrong variable in `PyramidNonrigidRegistration.cpp:54`"** | Verified: lines 54-55 store `downsampleFloatEnd` → `_downsampleFloatEnd` and `downsampleTargetEnd` → `_downsampleTargetEnd`. Code looks correct. Either already fixed or not-a-bug. **Closed — no action.** |
| **"Pattern A ignores `.flags` override"** | When a user passes a mesh object via Pattern A, `.flags` is read from the mesh. But there is no way to override flags independently. The `_prepare_arrays` extraction (above) includes `floating_flags` and `target_flags` override parameters. **Addressed by `_prepare_arrays` extraction.** |
| **"`RigidTransform.matrix` read-only"** | Currently `def_prop_ro` in `bindings.cpp`. Making it writable allows construction from user-supplied matrices. **Deferred to v0.3** — not blocking any v0.2 use case. |

## Branch operations checklist

The `meshmonk-modernization` branch is ~119 commits ahead of `origin/cli`. Execution order:

1. All v0.2 code/doc changes land on `meshmonk-modernization`
2. PR `meshmonk-modernization` → `cli` (merge commit — preserves history for the large change set)
3. PR `cli` → `master` (merge commit)
4. Rename `master` → `main` via GitHub Settings → Default Branch
5. Update any CI workflow branch filters (current `ci.yml` already triggers on both `["master", "main"]`)
6. Verify CI passes on `main`

README and migration guide must land **before** step 2, since the PR to `master` is the public-facing merge and visitors should see the correct README.

## README rewrite scope

Current README (`README.md`, 189 lines) has substantial stale content. The rewrite should:

**Remove:**
- C++ CLI instructions (lines 100-162) — ADR D4 deleted the C++ CLI
- Legacy `Makefile` build process
- `TheWebMonks/meshmonk` clone URL — repo transferred to `jsnyde0/meshmonk` (ADR D9)
- C++14 compiler requirements section

**Keep/update:**
- CI matrix table (update with current compilers)
- MATLAB redirect banner (verify fork URL)

**Add:**
- Installation via `pip install meshmonk` (or build-from-source instructions until PyPI)
- Python API quickstart (3-4 line example for rigid, nonrigid, pyramid)
- Python CLI usage summary
- Link to migration guide (`docs/migration-from-matlab.md`)
- Link to ADR for architectural decisions
- Development setup instructions (for contributors)

## Type stubs

- Auto-generate `_meshmonk_core.pyi` via `nanobind.stubgen` at build time (ADR D2 mentions nanobind's "first-class `stubgen`")
- Add `py.typed` marker file at `meshmonk/py.typed`
- Python-level wrappers in `__init__.py` already use standard typing (dataclasses, type annotations) — no additional stubs needed
- CI verification: add a mypy or pyright check on a minimal import example

## Windows CI

Replace the vague "if feasible" with a concrete assessment:

**Build dependencies (all Windows-compatible):**
- OpenMesh 11.0.0: vendored, CMake, supports MSVC
- Eigen 3.4.0: header-only
- nanobind: supports MSVC 19.30+ (VS 2022)
- scikit-build-core: supports Windows
- C++20 designated initializers: MSVC 19.30+ supports these

**Implementation checklist:**
1. Add `windows-latest` + MSVC 19.30 to CI matrix
2. If OpenMesh emits warnings-as-errors under MSVC, compile OpenMesh targets at C++17 (mitigation from overview design)
3. Success criterion: full test suite passes on Windows

## cibuildwheel dry-run

**Definition:** Add a CI job that builds one manylinux wheel using cibuildwheel. No PyPI upload.

**Success criteria:**
- `cibuildwheel --print-build-identifiers` shows expected targets
- One manylinux2014_x86_64 wheel builds successfully
- Built wheel passes `pytest` inside the cibuildwheel container

Full wheel matrix (all OSes, all Python versions) is v0.3 scope.
