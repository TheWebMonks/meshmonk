# Design: MeshMonk Modernization — Python-First Rebuild

**Date:** 2026-04-17
**Status:** Draft (ready for v0.0 implementation)
**Decisions:** [ADR-001](../docs/decisions/ADR-001-meshmonk-modernization.md)
**Origin:** Brainstorming session 2026-04-17; prompted by revisiting the May 2025 "Phase 5" spec and finding it too narrow
**Supersedes:** [`docs/specs/2025-05-29-modernization-and-python-bindings.md`](../docs/specs/2025-05-29-modernization-and-python-bindings.md)

---

## Problem

MeshMonk is a C++ 3D mesh registration library (rigid ICP, nonrigid deformable registration, pyramid multi-resolution). Algorithmic core represents ~10 years of research-validated tuning by a university lab. A 2025 modernization effort (Phases 1–4) reorganized the repo and added first-pass pybind11 bindings written by a then-young AI agent.

The current state has real structural problems that go beyond what the May 2025 "Phase 5" plan (demo scripts + README polish) would fix:

- **No packaging:** users `sys.path.append('./build/python')`; no `pip install`, no `pyproject.toml`, no wheels
- **Untrusted bindings:** the pybind11 wrapper was never validated numerically; the author's own comments flag confusion about the API shape
- **API smells:** 15-arg telescoping signatures, in-place mutation of numpy arrays, no typed failure mode, parallel `_mex` raw-pointer API surface
- **No test coverage:** the one test file (`test_bindings.py`) is a `print`-based smoke with random data; no CI; no numerical assertions. No C++ unit tests exist either (OpenMesh vendors its own, but MeshMonk core has none)
- **MATLAB burden:** ~45% of non-algorithmic code is MATLAB MEX plumbing, blocking cleanup
- **C++ CLI of unknown quality:** built by the same untrusted AI-agent run with no tests; pending verification work on branch `add-compute-rigid-transform-cli`

Meanwhile, parallel work streams (Karlijne's CLI verification, the MATLAB-CMake adaptation, and the Python bindings themselves) were all branching off `cli` without a coherent strategy for where the project is actually going.

## Goal

A **clean, modern, Python-first mesh registration library** the owner is proud of — not a patched-up version of the current state. Python is the primary interface for end users and AI agents; a minimal Python CLI replaces the C++ CLI; MATLAB support moves to a university fork; the C++ core is modernized to C++20 with a redesigned public API; every piece has a real harness so we iterate against reality, not interpretation.

## Non-Goals

Out of scope for v0.1 (the initial overhaul), explicit:

- Cross-platform wheels on PyPI — deferred to v0.3
- MCP server wrapping the library for agents — deferred to v0.4+
- Dropping OpenMesh in favor of `meshoptimizer` — deferred to v0.4+
- Meson migration — contingent on OpenMesh drop, v0.4+ at earliest
- Windows CI — likely v0.2, certain by v0.3
- Sphinx / MkDocs documentation site — v0.3
- GPU/CUDA backend — no timeline
- Algorithmic improvements to registration math — never, unless a measurable regression demands it
- Backward compatibility with the current `meshmonk_python` module name, `sys.path` hack, or pybind11-style API
- Salvaging the existing `test_bindings.py`, `example.cpp`, or the bulk of `add-compute-rigid-transform-cli` (verified 100% plumbing, zero library value). Exception: the 8-vertex cube fixture under `cli/test_data/rigid_transform/` is cherry-picked into `tests/fixtures/` before that branch is retired.

---

## Proposed Architecture

### Overview

```
┌─────────────────────────────────────────────────────────────┐
│  User (Python scripts, agents, CLI)                          │
└──────────────┬──────────────────────────────────────────────┘
               │  import meshmonk     $ meshmonk rigid in.obj target.obj
               ▼
┌─────────────────────────────────────────────────────────────┐
│  meshmonk/ (Python package, shipped on PyPI in v0.3)         │
│  ├── __init__.py     kwargs API, duck-typed mesh objects     │
│  ├── cli.py          typer CLI, entry point `meshmonk`       │
│  ├── io.py           optional trimesh/meshio/open3d helpers  │
│  ├── py.typed        PEP 561 marker                          │
│  └── _core.*.so      nanobind extension (compiled)           │
└──────────────┬──────────────────────────────────────────────┘
               │ nanobind (zero-copy read-only inputs; floating buffer copied once for mutation)
               ▼
┌─────────────────────────────────────────────────────────────┐
│  libmeshmonk  (C++20 shared library)                         │
│  Public API:                                                  │
│    rigid_registration(...)      → std::expected<Result,Err>  │
│    nonrigid_registration(...)   → std::expected<Result,Err>  │
│    pyramid_registration(...)    → std::expected<Result,Err>  │
│    compute_correspondences(...) → Correspondences            │
│    compute_inlier_weights(...)  → VecDynFloat                │
│    compute_rigid_transform(...) → std::expected<T,Err>       │
│    ...                                                        │
│  Deps: Eigen3, OpenMesh (→ meshoptimizer v0.4+), nanoflann   │
└─────────────────────────────────────────────────────────────┘
```

### Repository layout (after v0.1)

```
meshmonk/
├── pyproject.toml              # scikit-build-core + nanobind
├── CMakeLists.txt              # root: library + python
├── README.md                   # Python-first + MATLAB-fork redirect
├── library/                    # C++ library
│   ├── CMakeLists.txt
│   ├── include/meshmonk/
│   │   ├── meshmonk.hpp        # public API surface
│   │   ├── types.hpp           # canonical aliases (kill typedef duplication)
│   │   ├── params.hpp          # Rigid/Nonrigid/PyramidParams structs
│   │   ├── result.hpp          # Rigid/Nonrigid/PyramidResult structs
│   │   └── transform.hpp       # RigidTransform strong type
│   └── src/                    # internal implementation
├── python/
│   ├── src/
│   │   ├── _core.cpp           # nanobind bindings
│   │   └── meshmonk/           # Python package
│   │       ├── __init__.py
│   │       ├── cli.py
│   │       ├── io.py
│   │       └── py.typed
│   └── CMakeLists.txt
├── tests/
│   ├── conftest.py
│   ├── utils/
│   │   └── mesh_compare.py     # ~20-line numpy Hausdorff/RMSE helper
│   ├── fixtures/
│   │   ├── cube_rigid/         # ported from cli/test_data/rigid_transform
│   │   └── synthetic/          # generated ~200-vertex sphere pairs
│   ├── golden/                 # v0.1 human-approved goldens
│   ├── test_rigid.py
│   ├── test_nonrigid.py
│   ├── test_pyramid.py
│   ├── test_primitives.py
│   └── test_cli.py
├── data/                       # Template.obj, demoFace.obj (EXCLUDED from wheel; fetched on-demand)
├── vendor/                     # OpenMesh-11.0.0, nanoflann
├── docs/
│   └── decisions/              # ADRs
├── history/                    # design docs (this one, future ones)
└── .github/workflows/
    └── ci.yml
```

**Deleted from current state:** `cli/`, `demo/matlab/`, all `_mex` C++ functions, `test_meshmonk_mexing*`, old `test_bindings.py`, `library/examples/` (dead scratchpad with hardcoded absolute paths — Python CLI is the canonical reference; C++ consumers link against the installed library from their own CMake). Branches `feat-python-bindings`, `add-compute-rigid-transform-cli`, `feature/*`, `feat/meshmonk-cli*` retired.

### C++ API shape (v0.1)

High-level pipelines plus low-level composable primitives, all sharing the same param structs:

```cpp
namespace meshmonk {

    // Canonical types (library/include/meshmonk/types.hpp)
    using FeatureMat  = Eigen::Matrix<float, Eigen::Dynamic, 6>;
    using FacesMat    = Eigen::Matrix<int,   Eigen::Dynamic, 3>;  // signed int32 for numpy/trimesh interop
    using VecDynFloat = Eigen::VectorXf;
    inline constexpr int NUM_FEATURES = 6;

    // Strong rigid transform (library/include/meshmonk/transform.hpp)
    struct RigidTransform {
        Eigen::Matrix4f matrix = Eigen::Matrix4f::Identity();   // 4x4 SE(3); similarity if params.use_scaling
        RigidTransform compose(const RigidTransform&) const;    // returns (*this) * rhs
        FeatureMat     apply(const FeatureMat&) const;          // transforms positions; rotates normals
        RigidTransform inverse() const;
    };

    // Composable parameter structs (library/include/meshmonk/params.hpp)
    struct CorrespondenceParams { bool symmetric=true; int num_neighbours=3;
                                   float flag_threshold=0.9f; bool equalize_push_pull=false; };
    struct InlierParams { float kappa=12.0f; bool use_orientation=true; };
    struct ViscoElasticParams {
        float sigma = 3.0f;
        // Viscous + elastic annealing pairs interpolate across num_iterations.
        // Defaults mirror current MATLAB demos (start == num_iterations, end == 1).
        int num_viscous_iterations_start = 200;
        int num_viscous_iterations_end   = 1;
        int num_elastic_iterations_start = 200;
        int num_elastic_iterations_end   = 1;
    };
    struct DownsampleSchedule { float float_start=50.f, target_start=70.f,
                                      float_end=0.f,   target_end=0.f; };

    struct RigidParams    { CorrespondenceParams correspondences; InlierParams inliers;
                            int num_iterations=80; bool use_scaling=false; };
    struct NonrigidParams { CorrespondenceParams correspondences; InlierParams inliers;
                            ViscoElasticParams transform; int num_iterations=200; };
    struct PyramidParams  { /* + DownsampleSchedule; num_pyramid_layers=3; */ };

    // Result structs (library/include/meshmonk/result.hpp)
    // v0.1 ships minimal shapes; diagnostic fields (fitness, inlier_rmse,
    // num_inliers, converged) are deferred to v0.2 once detection sites exist.
    using Vec3Mat = Eigen::Matrix<float, Eigen::Dynamic, 3>;  // (N,3) positions or displacements

    enum class RegistrationError {
        DegenerateInput, DecompositionFailed, NonConvergence, InsufficientInliers
    };
    struct RigidResult {
        FeatureMat     aligned_features;   // (N,6) positions + normals
        RigidTransform transform;          // 4x4 SE(3)
        int            iterations_run = 0;
    };
    struct NonrigidResult {
        FeatureMat    aligned_features;      // (N,6)
        VecDynFloat   final_inlier_weights;  // (N,) weight per floating vertex
        Vec3Mat       displacement_field;    // (N,3) per-vertex displacement from original floating
        int           iterations_run = 0;
    };
    struct PyramidResult {
        FeatureMat         aligned_features;    // (N,6) at finest layer
        VecDynFloat        final_inlier_weights;
        Vec3Mat            displacement_field;  // (N,3) at finest layer
        std::vector<int>   per_layer_iterations;
    };

    // Public API (library/include/meshmonk/meshmonk.hpp)
    [[nodiscard]] std::expected<RigidResult, RegistrationError>
    rigid_registration(
        Eigen::Ref<const FeatureMat>  floating_features,
        Eigen::Ref<const FeatureMat>  target_features,
        Eigen::Ref<const FacesMat>    floating_faces,
        Eigen::Ref<const FacesMat>    target_faces,
        Eigen::Ref<const VecDynFloat> floating_flags,
        Eigen::Ref<const VecDynFloat> target_flags,
        const RigidParams&            params = {});

    // Same shape for nonrigid_registration, pyramid_registration.
    // Low-level primitives (compute_correspondences, compute_inlier_weights,
    // compute_rigid_transform, compute_nonrigid_transform, downsample_mesh,
    // scale_shift_mesh, compute_normals) share the same param structs.
    // (`compute_rigid_transform` / `compute_nonrigid_transform` are renames of
    // the old `compute_rigid_transformation` / `compute_nonrigid_transformation`.)
}
```

**Param struct field-order stability rule:** field order inside each param struct is STABLE across minor versions (v0.x → v0.x+1). Reordering fields requires an ADR update. The Python shim uses named `nb::arg()` bindings (not positional construction) so Python callers are immune to field reordering, but C++ callers using designated initializers would be. This protects against silent field-swap hazards (e.g., `num_iterations` and `num_neighbours` both being `int` and adjacent).

**Error detection sites (v0.1):** `RegistrationError` is a real enum, not aspirational. Current code writes `std::cerr` on failure and continues with garbage; the rewrite converts these into `std::expected` returns at four call sites:

| Error | Detected at | Criterion |
|---|---|---|
| `DegenerateInput` | boundary pre-check | empty inputs, dim/shape mismatch, `FeatureMat` rows ≠ `floating_flags.size()`, zero-range bbox, or all-zero inlier flags |
| `InsufficientInliers` | after `InlierDetector` (boundary-converted) | <4 non-zero inlier weights (rigid needs ≥3 for SVD well-posedness) |
| `DecompositionFailed` | internal (boundary-converted) | `RigidTransformer` SVD or `EigenVectorDecomposer` fails to converge on valid-shaped input |
| `NonConvergence` | wrapper loop (v0.2+ only) | convergence criterion (TBD) not met within `num_iterations`; **not raised in v0.1** since convergence tracking is deferred |

**Split between boundary and internal detection:** `DegenerateInput` is raised by an explicit pre-check at the top of each public API function (`rigid_registration()`, `nonrigid_registration()`, `pyramid_registration()`) — so shape / rank / dim errors never reach the algorithmic core. Internal algorithmic classes (`ViscoElasticTransformer`, `InlierDetector`, `RigidTransformer`) are **NOT** rewritten to return `std::expected`; they retain their existing control flow and signal failure via the log sink plus a thread-local status slot that the public wrapper reads and converts to `expected` before returning. This preserves the "load-bearing code untouched" constraint while still giving callers typed errors.

All internal `std::cerr` writes are routed through a logger sink (configurable via `meshmonk.set_log_level()` in Python) so consumers can suppress stderr noise. The log sink is thread-safe (atomic level read) but its effect is global to the process.

Load-bearing algorithmic code (`ViscoElasticTransformer`, `InlierDetector`, `RigidTransformer`, `NonrigidRegistration` annealing, `PyramidNonrigidRegistration` scheduling) remains untouched. The redesign is entirely at the API surface above them, plus confining OpenMesh to an I/O boundary (stop reconstructing `TriMesh` inside `update()` calls).

### Python API shape (v0.1)

```python
import meshmonk
import numpy as np
import trimesh  # optional

# 1. Pythonic kwargs + duck-typed mesh objects
mesh = trimesh.load("Template.obj")
result = meshmonk.rigid_register(
    floating=mesh,                 # any object with .vertices, .faces; also accepts path str
    target="demoFace.obj",
    num_iterations=80,             # kwargs override RigidParams defaults
    inlier_kappa=12.0,
)

result.aligned_vertices   # (N, 3) np.ndarray — positions only
result.aligned_features   # (N, 6) np.ndarray — positions + recomputed normals
result.transform.matrix   # (4, 4) np.ndarray
result.iterations_run     # int
# diagnostic fields (fitness, converged, inlier_rmse, num_inliers) → v0.2

# 2. Raw-arrays path (zero wrapping)
result = meshmonk.rigid_register(
    floating_features=V_float,     # (N,6) float32
    target_features=V_target,
    floating_faces=F_float,        # (M,3) int32
    target_faces=F_target,
    floating_flags=flags_float,    # (N,) float32
    target_flags=flags_target,
)

# 3. Low-level composable primitives
corrs = meshmonk.compute_correspondences(V_float, V_target, flags_float, flags_target)
weights = meshmonk.compute_inlier_weights(V_float, corrs, kappa=12.0)
transform = meshmonk.compute_rigid_transform(V_float, corrs, weights)

# 4. CLI (installed as `meshmonk` by pip entry point)
# $ meshmonk rigid    Template.obj demoFace.obj --out result.obj
# $ meshmonk nonrigid Template.obj demoFace.obj --out result.obj --iterations 200
# $ meshmonk pyramid  Template.obj demoFace.obj --out result.obj --layers 3

# 5. First-run experience (canonical smoke test)
# $ meshmonk demo --download          # fetches Template.obj + demoFace.obj → ~/.cache/meshmonk/
# $ meshmonk demo rigid               # runs registration, opens MeshLab/Blender if available,
#                                     # else writes result OBJ and prints the path
```

**Kwarg flattening rule:** every param struct field maps to a top-level kwarg using the prefix `<struct>_<field>`:

- `CorrespondenceParams` → `correspondences_symmetric`, `correspondences_num_neighbours`, `correspondences_flag_threshold`, `correspondences_equalize_push_pull`
- `InlierParams` → `inlier_kappa`, `inlier_use_orientation`
- `ViscoElasticParams` → `transform_sigma`, `transform_num_viscous_iterations_start`, etc.
- Fields on the top-level params struct (`RigidParams.num_iterations`, `RigidParams.use_scaling`, `NonrigidParams.num_iterations`, `PyramidParams.num_pyramid_layers`) stay **unprefixed** — they're the ergonomic top-level knobs.

Adding two fields with the same prefixed name across different sub-structs is a breaking change (requires an ADR update). Callers who prefer explicit struct-passing can still use `rigid_params=meshmonk.RigidParams(num_iterations=80)`.

**nanobind zero-copy caveats:** the zero-copy path requires `float32` for `FeatureMat` and `int32` for `FacesMat`, plus C-contiguous memory. Default `trimesh` arrays use `int64` faces — the wrapper does a one-time cast at the Python boundary so end users aren't penalized. numpy row-major layout is accepted into Eigen's default column-major storage via a stride-aware `Eigen::Ref<…, 0, Eigen::Stride<Dynamic,Dynamic>>` signature (≈1–2% overhead on KNN loops per our read of Eigen). If profiling flags this as hot in v0.2, internal storage switches to row-major — a larger refactor deferred.

**Normals handling:** precedence is (a) explicit `normals=V_normals` kwarg if passed; else (b) `mesh.vertex_normals` if the mesh object exposes the attribute and its L2 norm per row is non-zero; else (c) recomputed via `meshmonk.compute_normals(V, F)`. Pass `compute_normals=True` to force recomputation even when `.vertex_normals` is available. For the raw-array path, if `floating_features` columns 3–5 are all-zero (as current MATLAB demos pass — `single(zeros(size(floatingPoints)))`), the wrapper auto-recomputes normals and emits a one-line warning via the log sink, matching current C++ behavior while surfacing the smell for the migration. `meshmonk.features_from_vertices(V, F)` is the recommended helper to avoid the warning.

Errors bubble as Python exceptions (nanobind translates `std::expected` unexpected values via a registered translator). Type stubs (`.pyi`) shipped with the wheel from v0.2 onward.

### Python package dependencies

Hard runtime deps (in `[project] dependencies`): `numpy >= 1.24`, `typer >= 0.9`. Build-time requires (`[build-system] requires`): `scikit-build-core >= 0.8`, `nanobind >= 2.0`, `cmake >= 3.26`.

Optional I/O extras (in `[project.optional-dependencies]`):

- `meshmonk[io]` — pulls `trimesh >= 4.0` (chosen over `meshio` for richer mesh-object model, built-in normals computation, and wider format coverage). Used by the wrapper's path-input handling and by `meshmonk demo`.
- `meshmonk[dev]` — pulls `pytest`, `pre-commit`, `ruff` for contributors.

The CLI (`meshmonk demo`, `meshmonk rigid`, etc.) depends on the `io` extra since it loads OBJ files. The `meshmonk` entry point emits a clear `pip install 'meshmonk[io]'` hint if `trimesh` is missing.

### Staged roadmap

**v0.0 — Reality capture (~1–2 days):** establish trust anchors before refactoring. Create `overhaul/v0.0-prep` branch from `cli`. Cherry-pick the 8-vertex cube fixture from `origin/add-compute-rigid-transform-cli:cli/test_data/rigid_transform/` to `tests/fixtures/cube_rigid/` before that branch is retired. Move `demo/Template.obj`, `demo/demoFace.obj`, `demo/rigid_transform.txt`, `demo/{Demo,Template}FaceLandmarks.csv` to `data/`. **Capture legacy baseline:** build the current `cli` branch C++ library + CLI, run it on `Template.obj ↔ demoFace.obj`, save `rigid_output.obj`, `pyramid_output.obj`, and `rigid_transform.txt` under `tests/golden/legacy_baseline/` (scientific-equivalence reference — see Tier 3.5). Write the mesh-compare helper. Set up pytest scaffolding. Commit this design doc + ADR. No library changes.

**v0.1 — Foundation (~1–2 focused weeks):** C++20 overhaul + nanobind bindings + Python package + typer CLI + scikit-build-core + GitHub Actions CI. `uv pip install .` (non-editable) works on the owner's Mac; editable installs tracked via [astral-sh/uv#14383](https://github.com/astral-sh/uv/issues/14383) — fallback is `pip install -e . --no-build-isolation`. All Tier 1–5 harness tests passing, including **owner runs one Blender/MeshLab visual sign-off session for rigid / nonrigid / pyramid outputs** to freeze Tier 3 goldens. No `meshmonk::legacy::` hedge — the current public API has no external consumers to preserve (MATLAB/C++-CLI/pybind11 all being dropped).

- **CI matrix (v0.1):** Ubuntu 22.04 + Ubuntu 24.04; macOS 14 (arm64) + macOS 13 (x86_64); Python 3.10, 3.11, 3.12, 3.13; CMake 3.26+ installed via `lukka/get-cmake@latest` at a pinned version; gcc 11+, clang 15+, AppleClang 15+. Matrix is explicit in `.github/workflows/ci.yml` and mirrored in a CI-matrix table in the README.
- **Developer tooling (v0.1):** ships `.pre-commit-config.yaml` (ruff for Python lint+format, `clang-format-16` with LLVM style for C++, `cmake-format` for CMakeLists, `check-yaml`, `end-of-file-fixer`), `.editorconfig` (4-space Python, 4-space C++, LF line endings), and `CONTRIBUTING.md` stub. CI runs `pre-commit run --all-files` as a first-check job.

**v0.2 — Team-ready:** ship type stubs. README rewrite. Migration guide from MATLAB (side-by-side snippets) at `docs/migration-from-matlab.md`. Windows CI if feasible. cibuildwheel dry-run. Add diagnostic fields (`converged`, `fitness`, `inlier_rmse`, `num_inliers`) to result structs once convergence criteria are pinned down. PR `cli` → `master`; rename `master` → `main`.

**v0.3 — PyPI-ready:** cibuildwheel matrix (manylinux, macOS x86_64 + arm64, Windows). Test PyPI round-trip. MkDocs-material docs site. Publish to PyPI. `pip install meshmonk` works globally.

**v0.4+ — Future:** OpenMesh → meshoptimizer; reconsider Meson; MCP server; benchmarks.

### Branch strategy

```
master (2021, stale)
  ↑ PR at v0.2 milestone, then rename master → main
  │
cli (integration trunk — current)
  ↑ merge v0.1 overhaul PR
  │
overhaul/v0.0-prep  ← branched from cli, where v0.0 work happens
overhaul/v0.1       ← branched from v0.0-prep, where main rebuild happens
```

Branches to retire: all remote branches except `master`, `cli`, and the new `overhaul/*` — explicitly including `feat-python-bindings`, `add-compute-rigid-transform-cli` (after cube-fixture cherry-pick), `feature/matlab-cmake-adaptation`, `feature/phase1-cmake-restructure`, `feat/meshmonk-cli`, `feat/meshmonk-cli-superbuild-attempt`, `cpp_project`, `development`, `mh-fixdocs`.

### Canonical parameter defaults

Use the MATLAB-demo-tuned values as the Python API defaults (not the C++ header defaults, which were a scratchpad developer's guesses). Full table in [ADR-001, D6 Rationale](../docs/decisions/ADR-001-meshmonk-modernization.md#d6-redesign-the-public-api). Key differences from current C++ defaults:

| Param | Current C++ | MATLAB (adopted) |
|---|---|---|
| `inlier_kappa` | 4.0 | **12.0** |
| `correspondences_num_neighbours` | 5 | **3** |
| `correspondences_flag_threshold` (rigid/nonrigid) | 0.99 | **0.9** |
| `correspondences_flag_threshold` (pyramid) | 0.99 | **0.999** |
| `pyramid.num_iterations` | 60 | **90** |

### Convention: floating vs. target

**Enforced across entire API:** `floating` is the mesh being transformed (conceptually: the template we deform); `target` is the fixed reference (the face scan we're trying to match). Verified consistent across all MATLAB demos.

### Harness strategy

Five tiers, ordered fast → slow:

- **Tier 1 — Analytical** (deterministic, math-grounded, <5 s total): synthetic rigid recovery (apply known SE(3) to a mesh, run `rigid_registration` from identity, verify recovery), self-consistency, round-trip, correspondence sanity, inlier-weight sanity. This is where end-to-end ICP gets covered.
- **Tier 2 — Primitive fixture** (deterministic, <1 s): 8-vertex cube + known correspondences + inlier weights. Exercises `compute_rigid_transform` as a least-squares primitive, **not** end-to-end ICP. Ported from `cli/test_data/rigid_transform/`.
- **Tier 3 — Human-approved visual golden** (captured once per registration type in v0.1): owner inspects output in Blender or MeshLab, approves, we freeze as `.npz`. Future runs compare with Hausdorff/RMSE via `tests/utils/mesh_compare.py` (~20 lines of numpy).
- **Tier 3.5 — Legacy scientific-equivalence gate**: reference outputs captured in v0.0 on `Template.obj ↔ demoFace.obj` become the "did the rewrite preserve the current library's numerical behavior?" baseline. New-code output must match within calibrated tolerance. Trust comes from the *library* (research-validated per ADR D1), not from whatever plumbing wraps it. **Capture source, in priority order:** (a) the `cli`-branch C++ CLI binary if it builds cleanly on the owner's toolchain; (b) the MATLAB demos (`test_rigid_registration.m`, `test_pyramid_registration.m`) which drive the same research-validated library via MEX — equally valid as a baseline; (c) if neither is feasible in v0.0, Tier 3.5 becomes optional and the rewrite leans harder on Tier 1 + Tier 3 human sign-off. Not circular with ADR D7's rejection of untrusted pybind11 output — the baseline is either the raw library (a) or MEX into the raw library (b), not the abandoned binding layer.
- **Tier 4 — Soft anchor**: `data/rigid_transform.txt` (committed real-face rigid transform) — loose-tolerance sanity
- **Tier 5 — E2E smoke**: CLI commands exit 0, output file exists and is valid OBJ

**Tolerance calibration (v0.1 task):** Measure run-to-run noise across platforms (owner's Mac + Linux CI + macOS CI), run each registration 10× per platform, set Tier 3 / 3.5 tolerance at 2× the measured cross-platform max RMSE. Apply a **floor of 1e-4 mm absolute** — if the cross-platform max RMSE is below ~1e-5 mm (i.e., the code is deterministic within double rounding), the floor prevents tolerance collapse to zero and accommodates future compiler / libm / platform drift.

**Golden file format:** all Tier 3 and Tier 3.5 goldens are stored as `.npz` (float32 vertices, int32 faces, `(4,4)` float32 transform). A `tests/utils/obj_to_npz.py` converter ingests the legacy `.obj` captures during v0.0. Raw OBJs would inflate git (10+ MB per nonrigid golden); npz compresses to ~1.3 MB per 54K-vertex golden, keeping v0.1 comfortably under the threshold where LFS becomes necessary.

**Tier 3 runtime budget:** empirical estimate 30–180 s per registration type on `Template.obj` (7K verts) ↔ `demoFace.obj` (54K verts). Measure during v0.1; if the combined Tier 3 + 3.5 suite exceeds ~5 min per CI run per OS, move both to a nightly workflow and keep PR CI on Tiers 1, 2, 4, 5.

**Golden update protocol:** when an algorithm change legitimately requires updating a Tier 3 or 3.5 golden:
1. PR modifies algorithm AND updates `tests/golden/<name>.npz`
2. PR description includes before/after screenshot + RMSE delta against the old golden
3. Owner approves visually (reviewer sign-off not sufficient on its own)
4. Commit message includes `[GOLDEN-UPDATE]` tag
5. Old goldens are preserved via git history — no LFS needed at v0.1 scale (~MB-level npz)

**Visualization tooling:** outsourced to third-party (Blender, MeshLab). No in-repo viewer. Only the numpy compare helper is in-repo.

---

## Key Design Decisions

All decisions captured in [ADR-001](../docs/decisions/ADR-001-meshmonk-modernization.md):

- **D1 (FIRM):** Modernize in C++20 — not Rust
- **D2 (FIRM):** nanobind for Python bindings — not pybind11
- **D3 (FIRM):** Drop MATLAB support — redirect to university fork
- **D4 (FIRM):** Delete C++ CLI — replace with Python typer CLI
- **D5 (FIRM):** CMake + scikit-build-core — not Meson
- **D6 (FIRM):** Params structs + result structs + `std::expected`
- **D7 (FLEXIBLE):** 5-tier harness, analytical first, no circular baselining
- **D8 (FLEXIBLE):** Keep OpenMesh for v0.1; v0.4+ migration requires a half-edge-lite module, not just meshoptimizer
- **D9 (FIRM):** Transfer repo to `jsnyde0/meshmonk`; claim PyPI `meshmonk`

---

## Consequences

**What becomes easier:**
- Installing: `uv pip install meshmonk` (v0.3) vs. clone+build+sys.path hack
- Testing: CI-backed pytest covers library + CLI vs. no automation
- Writing agent-callable code: typed kwargs + result objects + exceptions vs. in-place-mutation guesswork
- Contributing: single API surface vs. three parallel ones (Eigen/`_mex`/pybind11)
- Extending the API: add a field to a struct vs. edit 15-arg signatures in four places

**What becomes harder:**
- MATLAB users must migrate to university fork (one-time cost for them)
- Contributors must know C++20 (minor barrier for new contributors; major improvement for correctness)
- Branch-cleanup dance: a few "retire" PRs needed before the repo history becomes clean

**Tradeoffs:**
- Deleting untrusted code and rebuilding costs calendar time vs. patching (~1–2 weeks) but produces something defensible
- Keeping OpenMesh defers a ~10–15% binary-size win; we accept that for v0.1 scope control
- Using `tl::expected` (header-only, single file, identical API to `std::expected`) unconditionally for v0.1 sidesteps the MSVC/libc++/libstdc++ `<expected>` version matrix; migrate to `std::expected` in v0.2 once the compiler matrix is known

**What systems need to change:**
- CMake: require `cmake >= 3.26` at the root (needed for `SKBUILD_*` variables and nanobind's `nanobind_add_module`). Replace `file(GLOB_RECURSE)` in `library/CMakeLists.txt` with explicit source enumeration (official CMake guidance — globs skip incremental rebuilds). `pyproject.toml` + `install(TARGETS … DESTINATION ${SKBUILD_PLATLIB_DIR}/meshmonk)` + delete CLI and MEX targets.
- Tests: from zero → pytest-based with fixtures, CI, and goldens. `tests/` lives at repo root, NOT shipped in the wheel — dev-only, invoked via `pytest tests/` from a repo checkout.
- Branch structure: retirement of 9+ obsolete branches plus two new `overhaul/*` branches
- README: completely rewritten, Python-first

### Release & versioning

- **v0.x** = pre-1.0; Python API may break at any minor bump
- **v1.0** = Python API freeze. C++ ABI explicitly NOT promised (Python-only users are stable; C++ direct consumers may need rebuilds)
- **Python & numpy floor:** v0.1 targets CPython 3.10–3.13. v0.3 wheels are built with nanobind's stable-ABI mode (`NB_STABLE_ABI=1`) against numpy 2.0 — forward-compatible with numpy 1.x at runtime. One wheel per OS/arch covers all supported Pythons (6 wheels total for v0.3: manylinux x86_64, manylinux aarch64, musllinux x86_64, macOS arm64, macOS x86_64, Windows x86_64 if included). If profiling shows abi3 overhead exceeds ~2%, fall back to per-Python wheels in v0.4 — no API impact, just a build-flag change.
- **Licensing continuity** (under D9 repo transfer): Apache-2.0 continues across the full tree. A `NOTICE` file is added at repo root attributing original copyright to TheWebMonks (2017) and new copyright to jsnyde0 (2026). `pyproject.toml` `[project]` lists `license = {text = "Apache-2.0"}`, `authors = [{name = "jsnyde0"}]`. New source files carry `Copyright 2026 jsnyde0` headers; untouched legacy files keep their original headers.
- **Build-time dependency strategy:** v0.1 vendors Eigen 3.4.0 under `vendor/eigen-3.4.0/` (~3 MB headers only) with SHA-256 pin in `vendor/VERSIONS.txt`. OpenMesh 11.0.0 remains vendored. nanoflann is already header-only-vendored. **No `FetchContent` at build time** — `pip install` works offline and cibuildwheel has no network dependency during compile.
- **sdist contents:** ships `library/`, `python/`, `vendor/`, `CMakeLists.txt`, `pyproject.toml`, `LICENSE`, `NOTICE`, `README.md`. Excludes `tests/`, `data/` (demo assets fetched separately), `history/`, `docs/`, `tutorial/`, `.github/`. Expected sdist size ~30 MB dominated by vendored OpenMesh + Eigen. `pip install meshmonk --no-binary :all:` works on any platform with a C++20 toolchain.
- **MATLAB migration**: README banner links to `docs/migration-from-matlab.md` (written in v0.2) — shows MATLAB demo snippets paired with equivalent Python kwargs calls
- **Release workflow** (from v0.3): GitHub Actions → cibuildwheel → PyPI via trusted-publisher; GitHub release tag and PyPI publish happen atomically in one workflow
- **PyPI trusted-publishing setup (v0.3 prerequisite):** before first publish: (a) claim the PyPI `meshmonk` name by uploading a placeholder `0.0.0.dev0` from `jsnyde0` manually; (b) configure trusted-publisher on `pypi.org/manage/account/publishing/` bound to `jsnyde0/meshmonk` + `.github/workflows/release.yml` + environment `pypi-release`; (c) dry-run on TestPyPI first. Workflow requires `permissions: id-token: write`.
- **Repository governance**: the repo transfers from `TheWebMonks/meshmonk` to `jsnyde0/meshmonk` before v0.1 lands (preserves the 85 stars, redirects existing `git clone` URLs, and anchors the Python rewrite under personal attribution). PyPI namespace `meshmonk` will be claimed by `jsnyde0`. If TheWebMonks continues MATLAB work, they fork from the transferred repo.
- **Documentation scope**: v0.1 = inline doxygen comments + README + ADRs; v0.3 = mkdocs-material site

---

## Known Limitations

**Explicitly not solving:**
- OpenMesh is still the halfedge/decimation dependency — binary size, build time, and transitive dep burden deferred to v0.4+. If OpenMesh 11.0.0 emits C++20 warnings-as-errors, the fallback is to compile the `OpenMeshCore` + `OpenMeshTools` subtargets at C++17 (via `target_compile_features(... cxx_std_17)`) while MeshMonk itself stays on C++20 — this works because OpenMesh types are internal-only, never in MeshMonk's public API surface.
- `tl::expected` (header-only) is used unconditionally for v0.1 to sidestep MSVC/libstdc++/libc++ `<expected>` version gaps; `std::expected` migration is scheduled for v0.2.
- No image-based visual regression tests — if a future nonrigid regression is visually subtle but not numerically detectable, Tier 3 may miss it. Mitigation: user does one-time visual sign-off; subsequent numerical comparison is strict, and the Tier 3.5 legacy-baseline gate catches most numerical drift.
- No GPU backend — MeshMonk remains CPU-only for v0.1 through v0.3
- **Windows compatibility is unverified end-to-end.** First Windows build attempt in v0.1 may surface MSVC `<expected>` gaps (mitigated by `tl::expected` above), OpenMesh 11 MSVC warnings, and OpenMP bundling issues. Fallback: ship manylinux + macOS wheels at v0.3 and defer Windows wheels to v0.4.

**Wheel distribution caveats:**
- **Demo data** (`Template.obj` ~820 KB, `demoFace.obj` ~10.6 MB) lives in `data/` for local development but is EXCLUDED from the wheel via `[tool.scikit-build.wheel.exclude]`. Users running `meshmonk demo` fetch it on-demand from a GitHub release asset, pinned by SHA-256, cached in `~/.cache/meshmonk/`.
- **OpenMP**: OpenMesh's CMake enables OpenMP by default, which transitively adds a `libomp` dependency. We either (a) disable OpenMP in the OpenMesh subbuild (`set(USE_OPENMP OFF)` before `add_subdirectory`), or (b) bundle `libomp` in the wheel via `delocate-wheel` (macOS) / `auditwheel repair` (Linux). Decision in v0.1.

**Concurrency gotchas:**
- The library is single-threaded for v0.1 (parallel correspondence search via OpenMP may be a later optimization). Python GIL release is easy to add in nanobind and should be included from v0.1 to avoid blocking the Python interpreter during long registrations.
- **Thread-safety contract:** public API functions are reentrant — concurrent calls on *disjoint* inputs from different Python threads are safe. The library holds no mutable module-level state except the log sink, which uses atomic level access. `meshmonk.set_log_level()` is thread-safe but its effect is global to the process. nanobind releases the GIL around each registration call so `concurrent.futures.ThreadPoolExecutor` scales across cores when the library is compiled with OpenMP (v0.2+).

---

## Open Questions

1. **Tolerance thresholds for Tier 3 / 3.5 goldens** — needs empirical measurement during v0.1. Plan: run each registration 10× with identical inputs across platforms (owner's Mac + Linux CI + macOS CI), measure max cross-platform RMSE, set tolerance at 2× measured max with a 1e-4 mm absolute floor.

2. **Cross-platform CI matrix start date for Windows** — Ubuntu + macOS in v0.1 is confirmed. Windows could be v0.1 (if straightforward) or v0.2 (if friction). Decide during v0.1 based on first Windows build attempt.

3. **MATLAB university-fork URL for README redirect** — the owner to provide the canonical fork URL (the third-party university team's copy) before v0.2.

4. **Canonical paper citation** — MeshMonk is academic work. Which paper does one cite? Needed for README, docs, and any downstream research use.

5. **Redistribution rights for bundled mesh data** — `data/demoFace.obj` header comment attributes it to "3Q Technologies Ltd." (`support@3q.com`, unreachable in 2026). License terms unknown. `tutorial/TutorialData/SimulatedMappedFaces/SIMPOP_*.obj` similar unknown origin. Resolve before v0.3 PyPI publication; if restricted, host externally via GitHub release asset (which is already the plan for wheel size reasons — see Wheel distribution caveats). If restrictive, regenerating Tier 3 goldens against a licensed-clean replacement triggers the Golden update protocol.

---

## Ready-to-execute v0.0 kickoff checklist

When execution starts:

- [ ] Create branch `overhaul/v0.0-prep` from `cli`
- [ ] Verify `cli` builds on macOS with current Xcode + CMake — document state only, don't fix regressions
- [ ] **Capture legacy baseline (Tier 3.5)**: preferred path — build the `cli`-branch C++ library + `meshmonk_cli`, run `rigid_reg` and `pyramid_reg` on `Template.obj ↔ demoFace.obj`, save outputs under `tests/golden/legacy_baseline/` as `.npz` (via `tests/utils/obj_to_npz.py`). **Fallback**: if `cli` doesn't build on current Xcode/CMake, run MATLAB demos (`test_rigid_registration.m`, `test_pyramid_registration.m`) instead — same research-validated library via MEX. Document which path was taken.
- [ ] Cherry-pick `cli/test_data/rigid_transform/` from `origin/add-compute-rigid-transform-cli` into `tests/fixtures/cube_rigid/` — **only these files**: `input_vertices.txt`, `input_mesh.obj`, `corresponding_features.txt`, `inlier_weights.txt`, `expected_vertices.txt`, `expected_transform.txt`. Exclude `output_mesh_generated.obj` and `output_transform_generated.txt` (runtime artifacts from the old CLI, not fixture inputs).
- [ ] Move `demo/Template.obj`, `demo/demoFace.obj`, `demo/rigid_transform.txt`, `demo/{Demo,Template}FaceLandmarks.csv` → `data/`
- [ ] **Clean up `tutorial/` (99 MB)**: (a) retire `tutorial/TutorialScripts/*.m` alongside `demo/*.m` and `matlab/` under D3 (tag `pre-modernization` before deletion); (b) move `tutorial/TutorialData/SimulatedMappedFaces/SIMPOP_*.obj` → `data/simulated_faces/` (drop the 45 `.mtl` files — unused by the algorithm); (c) delete `tutorial/TutorialData/Template_Rigidly_Aligned.*` and `tutorial/TutorialData/Targets/demoFace.*` (duplicates of demo/ assets); (d) defer `TutorialSlides.pptx` until a licensed-clean alternative or owner-held slides exist. Target: ≥80 MB removed from working tree.
- [ ] Delete `library/examples/example.cpp` and the `library/examples/` directory (Python CLI is the canonical reference; C++ consumers link against the installed library).
- [ ] Set up `tests/` root with pytest scaffolding (`conftest.py`, empty `__init__.py`, minimal `pyproject.toml` or `pytest.ini`)
- [ ] Write `tests/utils/mesh_compare.py` — Hausdorff / RMSE / per-vertex distance helpers (~20 lines, numpy only)
- [ ] Write `tests/utils/obj_to_npz.py` — OBJ → npz converter for Tier 3.5 legacy-baseline ingestion
- [ ] Write first Tier 1 analytical test as a failing placeholder (synthetic rigid recovery)
- [ ] Add `.pre-commit-config.yaml` (ruff, clang-format-16, cmake-format), `.editorconfig`, `NOTICE` file (WebMonks 2017 + jsnyde0 2026 dual-copyright)
- [ ] Commit this design doc and ADR-001
- [ ] Open PR for v0.0 skeleton — review, merge to `cli`
- [ ] **Repo transfer**: initiate GitHub repo transfer `TheWebMonks/meshmonk` → `jsnyde0/meshmonk` (owner-controlled; preserves stars, forks, issues; old URL auto-redirects via GitHub)

