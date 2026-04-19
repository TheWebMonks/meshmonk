# Fixes: meshmonk-modernization-v0.2
Date: 2026-04-19
Review passes: 2 (architecture + implementation reviewers per pass)

## Critical
(none)

## Important

- **meshmonk/__init__.py:672** — `rigid_params` silently accepts non-dict truthy values (e.g., `True`, `42`). Add `TypeError` validation: `if not isinstance(rigid_params, dict): raise TypeError(...)`. Same fix at line 760.

- **tests/test_silent_logger.py:18** — Hardcoded `"python3.11"` in `LD_LIBRARY_PATH` construction. Breaks on Python 3.10/3.12/3.13 CI runners. Fix: `f"python{sys.version_info.major}.{sys.version_info.minor}"`.

- **meshmonk/__init__.py:265** — Normal validity check `np.any(norms > 0)` is too permissive: a mesh with 1 valid normal and 53,999 zero normals passes. Fix: change to `np.all(norms > 0)`.

- **library/include/meshmonk/params.hpp:62-64** — `PyramidParams` user-defined constructor makes it non-aggregate, breaking C++20 designated initializer contract (ADR D6). Fix: remove constructor, set `correspondences.flag_threshold` default inline or override at Python binding level.

- **library/include/meshmonk/types.hpp:13 + global.hpp:6** — Dual `NUM_FEATURES` constants in different namespaces (`meshmonk::NUM_FEATURES` and `registration::NUM_FEATURES`) with no compile-time linkage. Fix: add `static_assert(meshmonk::NUM_FEATURES == registration::NUM_FEATURES)` in `meshmonk.cpp`.

## Minor

- **README.md:72** — Says "pybind11 Python bindings" but nanobind is used (ADR D2 FIRM conflict). Fix: change to "nanobind Python bindings."

- **README.md:59-66** — CI matrix table omits Windows. Fix: add row `| Windows latest | 3.12 | MSVC |`.

- **library/src/meshmonk.cpp:161-162,254,363** — Stale `TODO (v0.2)` comments for convergence criterion, which was explicitly deferred to v0.3. Fix: change to `TODO (v0.3)`.

- **bindings/bindings.cpp:50-63** — `unwrap_expected` switch has no `default` case. If a new `RegistrationError` variant is added, throws with empty message. Fix: add `default: msg = "Unknown RegistrationError"; break;`.

## ADR Updates

- **ADR-001 D6**: No update needed. The `PyramidParams` constructor issue is a code fix, not a decision revision. The designated-initializer intent in D6 remains correct.
- **ADR-001 D2**: No update needed. The README "pybind11" reference is a docs error, not a decision discrepancy.

## Discarded

Findings triaged out of the fix list (autonomous decisions):

- **MeshMonkError.code not accessible from Python** (found by both reviewers): nanobind's `nb::exception` doesn't support binding extra members directly. String parsing via `"DegenerateInput" in str(e)` works at v0.2. Non-trivial fix requiring custom exception translator. Deferred to v0.3.
- **Duplicated kwargs-application logic (~107 lines)**: Refactoring scope. No functional impact. Defer to v0.3.
- **Duplicated validation boilerplate (~69 lines)**: Same rationale. Defer to v0.3.
- **integration-test CI job blocked by Windows failure**: CI architecture improvement. Windows is new; teething issues expected. Defer.
- **_REPO_DATA_DIR / demo mesh bundling**: Demo --download already deferred to v0.3. Installed users can't run demo regardless.
- **CI stub drift detection**: Stub is committed and build-generated. Drift detection is a v0.3 CI improvement.
- **displacement_field semantics with rigid_params**: Documentation gap. The behavior (nonrigid-only displacement) is correct; just needs docstring clarification in v0.3.
- **RigidTransformer zero sumWeights division**: Legacy code; post-hoc NaN check catches it. ADR D6 constrains internal class rewrites.
- **Uninitialized pointer members**: Legacy code, no observed bug. Adding `= nullptr` is safe but low priority.
- **Redundant imports in _prepare_arrays**: Cosmetic, no runtime impact.
- **Misleading demo --download message**: Related to v0.3 demo work.
- **test_warning_level_default lacks output assertion**: Test quality; not blocking.
- **No [tool.scikit-build] config**: scikit-build-core defaults are correct.
- **log() writes without synchronization**: v0.3 pluggable logger will address.
- **iterations_run returns configured value**: Correct today; v0.3 convergence will revisit.
- **Migration guide vague parameter mappings**: Minor docs improvement.
- **Design doc specifies mypy/pyright CI, not implemented**: Note as deferred in design doc.
- **Windows CI single Python version**: Acceptable for initial Windows support.
