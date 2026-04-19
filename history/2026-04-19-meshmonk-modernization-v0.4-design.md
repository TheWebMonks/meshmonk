# Design: MeshMonk v0.4 — API Polish

**Date:** 2026-04-19
**Status:** Draft
**Parent:** [Overview](./2026-04-17-meshmonk-modernization-design.md)
**Decisions:** [ADR-001](../docs/decisions/ADR-001-meshmonk-modernization.md)

---

## Scope

v0.4 collects deferred API polish and code quality improvements that were explicitly
postponed during v0.2 and v0.3 to stay focused on functional milestones. None of
these items add new registration algorithms. They fall into three buckets:

- **Observability** — diagnostic fields in result structs (`converged`, `fitness`,
  `inlier_rmse`, `num_inliers`); the `rigid_registration` faces parameter revisit
- **API surface** — `RigidTransform.matrix` writable; pluggable logger sink
- **Code quality** — `std::expected` migration; validation boilerplate dedup;
  CI stub drift detection; target normals check

Most items are independent and can be done in parallel on separate worktrees.
The diagnostic fields item is the largest and requires a research/spec pass before
any code is written.

---

## Provenance

| Item | First deferred | Also deferred |
|------|---------------|---------------|
| Diagnostic fields | v0.2 design doc line 14 | v0.3 disposition table line 26 |
| `std::expected` migration | v0.2 design doc line 16 | v0.3 disposition table line 27 |
| Pluggable logger sink | v0.2 design doc line 17 | v0.3 disposition table line 28; v0.2 fixes file line 52 |
| `RigidTransform.matrix` writable | v0.1 fixes file line 65 | v0.3 disposition table line 29 |
| Validation boilerplate dedup | v0.2 fixes file line 41 | v0.3 disposition table line 30 |
| CI stub drift detection | v0.2 fixes file line 44 | v0.3 disposition table line 33 |
| Target normals check | v0.1 fixes file line 29 | — |
| `rigid_registration` faces revisit | v0.1 fixes file line 68 | — |

---

## Items

### 1. Diagnostic fields (`converged`, `fitness`, `inlier_rmse`, `num_inliers`)

**Status:** Needs research before implementation
**Deferred from:** v0.2 design doc (line 14), v0.3 disposition table (line 26)
**Dependencies:** Defines convergence semantics — may interact with the
`NonConvergence` error code stub in `RegistrationError` (currently reserved for
v0.2+); also interacts with item 8 (`rigid_registration` faces parameter, since
faces are documented as "needed for future convergence check").

The legacy C++ registration functions run a fixed number of iterations. There is no
existing convergence criterion and no codebase definition of "fitness." Adding the
fields with undefined or arbitrary semantics would be worse than omitting them — the
v0.2 design doc used exactly this reasoning for the deferral.

**Research questions that must be answered before implementation:**

1. **`converged`** — What constitutes convergence? Options:
   - Transform delta below a threshold between iterations (change in translation
     magnitude + rotation angle)?
   - Inlier-weighted residual below a threshold?
   - The legacy code ran fixed iterations and called it done; is a convergence
     criterion a semantic addition or a backward-compatibility risk?

2. **`fitness`** — No codebase definition exists. Likely candidates from ICP
   literature:
   - Fraction of floating points that are inliers (inlier ratio)
   - Mean inlier-weighted correspondence residual
   - Something else specific to MeshMonk's Mahalanobis weighting?
   The definition must be pinned to a formula before the field is shipped.

3. **`inlier_rmse`** — RMSE of inlier-weighted correspondences after final
   iteration. Relatively straightforward once "inlier" is defined (see below), but
   the weighting scheme (Mahalanobis vs. binary threshold) must be stated.

4. **`num_inliers`** — `InlierDetector` produces continuous weights in [0, 1],
   not a binary flag. What threshold on these weights constitutes an "inlier" for
   counting purposes? A hard threshold (e.g., weight > 0.5) is the simplest choice,
   but the right value is empirical and may need to be exposed as a parameter.

**Proposed approach:**

- Write a research note (can be a section of the implementation bead) that proposes
  concrete definitions for all four fields, referencing ICP literature and the
  existing `InlierDetector` / `ViscoElasticTransformer` internals.
- Get owner sign-off on definitions before writing any C++ or Python.
- Add `converged` flag and `NonConvergence` error code to `RegistrationError` at
  the same time (the error code is already reserved but not raised).
- Ship as a single PR covering C++ structs + nanobind exposure + Python type stubs
  + tests.

**Scope note on ADR D6:** ADR-001 D6 states "Shipping fields with undefined semantics
would be worse than omitting them." This item exists precisely to resolve that
prerequisite — define the semantics first, then implement.

---

### 2. `RigidTransform.matrix` writable

**Status:** Ready to implement
**Deferred from:** v0.1 fixes file (line 65), v0.3 disposition table (line 29)
**Dependencies:** None

Currently `RigidTransform.matrix` is read-only (getter only; no setter). Users who
want to supply their own 4x4 homogeneous transform matrix — e.g., to initialize
registration from an external alignment or to compose with a matrix from another
library — have no clean way to do so.

**Proposed approach:**

- Add a C++ constructor `RigidTransform(const Eigen::Matrix4f& m)` that validates
  the input (determinant close to 1.0, last row `[0,0,0,1]`) and stores it.
- Expose a Python property setter `RigidTransform.matrix = np.ndarray` (4x4 float32)
  with the same validation, raising `ValueError` on invalid input.
- Alternatively, a static factory `RigidTransform.from_matrix(m)` is less surprising
  than a setter on a property and aligns with the `from_*` convention common in
  Python libraries. Recommend `from_matrix` over a property setter.
- Add a Tier 1 round-trip test: construct from matrix, apply, recover.

**Open question:** Should invalid matrices raise immediately (eager validation) or
lazily on first use? Eager validation is recommended — fail fast.

---

### 3. Pluggable logger sink

**Status:** Ready to implement (design to be confirmed)
**Deferred from:** v0.2 design doc (line 17), v0.3 disposition table (line 28);
v0.2 fixes file line 52 notes "`log()` writes without synchronization: v0.3
pluggable logger will address."
**Dependencies:** None (but thread-safety note below)

`meshmonk::log()` currently writes to `std::cerr`. `set_log_level("silent")`
suppresses output, which covers the primary use case. However, some users want to
capture logs programmatically — route them to a Python logging handler, a file, or a
test buffer.

The v0.2 fixes note also flags that `log()` writes without synchronization. A
pluggable sink implementation must address this: the sink registration must be
thread-safe (atomic or mutex-guarded), and the documentation must state whether the
sink callback itself is called under a lock (blocking concurrent registrations) or
without (requires the user's callback to be reentrant). Recommendation: call without
lock; document the requirement.

**Proposed approach:**

- Define a C++ sink type: `using LogSink = std::function<void(LogLevel, std::string_view)>`.
- Add `meshmonk::set_log_sink(LogSink)` and `meshmonk::reset_log_sink()` (resets to
  `std::cerr`).
- Guard the global sink pointer with `std::atomic` or a mutex.
- Expose `meshmonk.set_log_sink(callable)` in Python. The callable receives
  `(level: str, message: str)`. GIL acquisition must be explicit in the nanobind
  callback wrapper (since registration calls release the GIL).
- A `None` argument to `set_log_sink` resets to default (`std::cerr`).
- Add a test that installs a list-appending sink and verifies messages arrive.

**Interaction with `set_log_level`:** `set_log_level("silent")` suppresses at the
call site (no call to the sink). This behavior is preserved: the level gate fires
before the sink is invoked.

---

### 4. `std::expected` migration

**Status:** Ready to implement (compiler floor decision required first)
**Deferred from:** v0.2 design doc (line 16), v0.3 disposition table (line 27)
**Dependencies:** Requires a deliberate CI matrix change; no functional dependency
on other items.

`tl::expected` (header-only, single-file, Apache-2.0) provides an identical API to
C++23 `std::expected` and currently works correctly on all supported compilers.
Migrating to `std::expected` is a code quality / standards-compliance improvement,
not a functional one. The migration is a global rename plus a header swap — no logic
changes.

**Compiler floor implications:**

| Compiler | `std::expected` available since | Notes |
|----------|--------------------------------|-------|
| gcc | 13 (C++23 mode) | gcc-11, gcc-12 lack it |
| clang | 16 (C++23 mode) | clang-15 lacks it |
| MSVC | 19.36 (VS 17.6) | Available on current VS |
| AppleClang | 15 | Available |

Migrating drops Ubuntu 22.04 support (ships gcc-11 by default; upgrading to gcc-13
via `ubuntu-toolchain-r/test` PPA is possible but adds CI complexity). The v0.1
design doc acknowledged this tradeoff explicitly.

**Decision needed:** Drop Ubuntu 22.04 from the CI matrix, or add the gcc-13 PPA.
Recommendation: drop Ubuntu 22.04 — it reaches EOL in April 2027 and its Python
3.10 coverage is already supplied by `ubuntu-latest` (24.04) + explicit Python
version pinning.

**Proposed approach:**

- Update `CMakeLists.txt`: set `CMAKE_CXX_STANDARD 23` (or add `std::expected`
  availability check via `cmake_check_source_compiles`).
- Replace `#include "tl/expected.hpp"` with `#include <expected>` and
  `tl::expected<T,E>` → `std::expected<T,E>` across all library headers and source.
- Remove `vendor/tl-expected/` subtree.
- Update CI matrix: remove Ubuntu 22.04 row; update minimum compiler notes.
- No Python-layer changes (nanobind shim works with either).

---

### 5. Validation boilerplate dedup (~69 lines)

**Status:** Ready to implement
**Deferred from:** v0.2 fixes file (line 41), v0.3 disposition table (line 30)
**Dependencies:** None

`rigid_registration`, `nonrigid_registration`, and `pyramid_registration` each
contain near-identical input validation blocks (~23 lines each, ~69 total). The
duplication means a bug fix or a new validation rule must be applied in three places.

**Proposed approach:**

- Extract a `validate_registration_inputs(...)` free function in an internal header
  (not part of the public API surface) returning
  `tl::expected<void, RegistrationError>` (or `std::expected` post-item 4).
- Parameters: `floating_features`, `floating_flags`, `floating_faces`,
  `target_features`, `target_flags`, `target_faces` (the full input set; each
  registration function passes what it has).
- Call from all three registration functions before the per-function logic.
- Regression risk is low: the extracted function has no side effects; tests for
  each registration function's error paths cover it transitively.
- Consider a short table in a comment documenting which checks are universal vs.
  per-function (e.g., the `InsufficientInliers` check after the first inlier pass
  is rigid-specific and stays in place).

---

### 6. CI stub drift detection

**Status:** Ready to implement
**Deferred from:** v0.2 fixes file (line 44), v0.3 disposition table (line 33)
**Dependencies:** None

The CI configuration (`.github/workflows/ci.yml`) and `pyproject.toml` can drift
out of sync — for example, the Python versions tested in CI may diverge from the
`requires-python` constraint or the classifiers in `pyproject.toml`. Currently this
drift is only caught by human review.

**Proposed approach:**

- Add a lightweight Python script (`scripts/check_ci_matrix.py`) that:
  - Parses the Python versions from the CI matrix YAML (under `strategy.matrix`).
  - Parses `requires-python` and `[project.classifiers]` from `pyproject.toml`.
  - Asserts consistency: every Python version in the CI matrix is within the
    `requires-python` range and has a matching `Programming Language :: Python :: X.Y`
    classifier.
  - Exits non-zero with a descriptive message on mismatch.
- Add a `pre-commit` hook and/or CI job that runs this script on every push.
- Scope: Python versions only for now; could extend to OS matrix vs. cibuildwheel
  platforms in a follow-up.

---

### 7. Target normals check

**Status:** Ready to implement
**Deferred from:** v0.1 fixes file (line 29)
**Dependencies:** None (but logically grouped with item 5, validation boilerplate dedup)

The input validation in Python (`meshmonk/__init__.py`) checks whether the
*floating* features have all-zero normals and auto-recomputes them if so. The same
check is NOT applied to the *target* features. A target mesh with all-zero normals
will produce silently wrong correspondence results.

From v0.1 fixes (line 29): "Target features all-zero normals not auto-recomputed in
Pattern B. Only floating is checked. Fix: apply the same all-zero check +
auto-recompute to `feat_target` using `target_faces`."

Additionally, the C++ validation layer checks `target_flags.maxCoeff() == 0.0f` for
the flags array but has no equivalent guard for `target_features` normals columns
being all-zero (the floating normals check is at the C++ layer; the target normals
check exists only in Python).

**Proposed approach:**

- In `meshmonk/__init__.py` (Pattern B path): apply the same all-zero normal
  detection + auto-recompute to `feat_target` as is already done for `feat_floating`.
- In `library/src/meshmonk.cpp` input validation: add an all-zero check on the
  target normals columns (columns 3–5 of `target_features`), consistent with the
  existing floating-normals check.
- Add a test: target mesh with zero normals, verify auto-recompute fires and result
  is not degenerate.
- This item can be merged into item 5 (validation boilerplate dedup) if they are
  implemented together — the shared validator can include the target normals check.

---

### 8. `rigid_registration` faces parameter revisit

**Status:** Needs research; blocked on item 1 (diagnostic fields)
**Deferred from:** v0.1 fixes file (line 68)
**Dependencies:** Item 1 (diagnostic fields / convergence)

From v0.1 fixes file (line 68): "`rigid_registration` validates `faces` but doesn't
use them — API consistency design choice. Faces needed for future convergence check."

Currently `rigid_registration` accepts `floating_faces` and `target_faces` for API
shape consistency with `nonrigid_registration` and `pyramid_registration`, but the
rigid algorithm does not use face topology. The docstring should clarify this, but
the deeper question is: once convergence checking lands (item 1), will the rigid
convergence check need face topology?

**Open questions:**

- Does any plausible rigid convergence criterion require face topology? (Surface
  normal recomputation mid-iteration would, for example.)
- If no convergence criterion ever needs faces in the rigid path, should the API
  drop `floating_faces`/`target_faces` from `rigid_registration`? (Breaking API
  change — requires an ADR update and a migration note.)
- If faces will be used, keep them and improve the docstring now.

**Recommendation:** Resolve this as part of the item 1 research note. Do not change
the API surface until the convergence criterion is defined. Add a docstring
clarification in the meantime (tiny, safe, non-breaking).

---

## Parallelism

Items in **Group A** are independent of each other and of Group B. They can each run
on a separate worktree simultaneously:

**Group A — fully independent:**

| Item | Estimated scope | Notes |
|------|----------------|-------|
| 2. `RigidTransform.matrix` writable | Small | C++ + nanobind + 1 test |
| 3. Pluggable logger sink | Medium | C++ + nanobind + GIL care + tests |
| 5. Validation boilerplate dedup | Small | Refactor; covered by existing tests |
| 6. CI stub drift detection | Small | New script + pre-commit hook |
| 7. Target normals check | Small | Can fold into item 5 |

**Group B — interacting items (do in order or together):**

| Item | Dependency |
|------|-----------|
| 1. Diagnostic fields | Needs research note first; defines `converged`/`NonConvergence` |
| 8. `rigid_registration` faces | Blocked on item 1 research note |

**Group C — compiler floor change (sequence-sensitive):**

| Item | Dependency |
|------|-----------|
| 4. `std::expected` migration | Logically independent, but if done while Group A items are in flight on other worktrees, each worktree picks up the header change as a merge. Recommend doing item 4 first or last to minimize rebase churn. |

---

## Compiler floor change

Item 4 (`std::expected` migration) raises the effective compiler floor to
**gcc-13 / clang-16 / AppleClang-15 / MSVC 19.36**. Concretely:

- **Ubuntu 22.04 dropped** from CI matrix (ships gcc-11; no `std::expected`).
  Ubuntu 22.04 reaches EOL April 2027. Python 3.10/3.11 coverage continues via
  Ubuntu 24.04 + explicit Python version pins.
- **CI matrix after migration:** Ubuntu 24.04 (gcc-13+), macOS arm64
  (AppleClang-15+), Windows (MSVC 19.36+).
- **`vendor/tl-expected/` removed** — ~600-line single header; minor sdist
  size reduction.
- **No functional change** — `std::expected<T,E>` and `tl::expected<T,E>` have
  identical API for the subset MeshMonk uses.

This change is independent of all other v0.4 items. It can be done on its own
worktree and merged before or after the Group A items. The only coordination needed:
every open worktree picks up a modified `CMakeLists.txt` and removed vendor header
on merge; rebuild required but no logic changes.

---

## Risks

1. **Diagnostic fields scope creep** — Item 1 is the only item with genuine
   unknowns. If the convergence criterion research reveals that a real convergence
   check requires significant algorithmic changes (e.g., adding early-exit logic
   to the iteration loop), that work must be scoped separately. The ADR D6 constraint
   ("internal algorithmic classes are NOT rewritten") applies: a convergence check
   at the public wrapper level (comparing iteration deltas) is within scope; rewriting
   `ViscoElasticTransformer` internals is not.

2. **`std::expected` build breaks on unexpected compilers** — The CI matrix covers
   the known-supported compilers, but contributors using older toolchains will see
   build failures after the migration. The CONTRIBUTING.md should document the new
   minimum compiler versions. Users of the wheel are unaffected (prebuilt binary).

3. **Pluggable logger GIL interactions** — The nanobind callback wrapper for a
   Python log sink must reacquire the GIL before calling the Python callable (since
   registration functions release the GIL). Forgetting this causes a crash, not a
   test failure. The implementation must include a test that exercises the sink from
   a `ThreadPoolExecutor` call (not just single-threaded) to catch GIL errors.

4. **API breaks from item 2 or 8** — `RigidTransform.from_matrix()` is additive
   (no break). If item 8 decides to drop `floating_faces`/`target_faces` from
   `rigid_registration`, that is a breaking API change requiring an ADR update and
   migration note. No break should happen without explicit owner sign-off.

5. **Item 5 regression risk is low but non-zero** — Extracting shared validation
   logic is a mechanical refactor, but any subtle behavioral difference (e.g., check
   order, short-circuit logic) could silently change which error code is returned for
   an edge-case input. The existing error-path tests provide a safety net; a
   validation-specific test suite covering all error cases for all three registration
   functions is recommended before the dedup PR lands.
