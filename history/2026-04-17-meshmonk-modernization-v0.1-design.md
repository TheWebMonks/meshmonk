# Design: MeshMonk v0.1 — Foundation

**Date:** 2026-04-17
**Status:** Draft
**Parent:** [Overview](./2026-04-17-meshmonk-modernization-design.md)
**Decisions:** [ADR-001](../docs/decisions/ADR-001-meshmonk-modernization.md)

---

## Scope

**v0.1 — Foundation (~1–2 focused weeks):** C++20 overhaul + nanobind bindings + Python package + typer CLI + scikit-build-core + GitHub Actions CI. `uv pip install .` (non-editable) works on the owner's Mac; editable installs tracked via [astral-sh/uv#14383](https://github.com/astral-sh/uv/issues/14383) — fallback is `uv pip install nanobind scikit-build-core cmake ninja && pip install -e . --no-build-isolation` (build-isolation off means build deps must be in the active env). All required harness tiers passing, including **owner runs one Blender/MeshLab visual sign-off session for rigid / nonrigid / pyramid outputs** to freeze Tier 3 goldens. No `meshmonk::legacy::` hedge — the current public API has no external consumers to preserve (MATLAB/C++-CLI/pybind11 all being dropped).

- **CI matrix (v0.1):** Ubuntu 22.04 + Ubuntu 24.04; macOS 14 (arm64) + macOS 13 (x86_64); Python 3.10, 3.11, 3.12, 3.13; CMake 3.26+ installed via `lukka/get-cmake@latest` at a pinned version; gcc 11+, clang 15+, AppleClang 15+. Matrix is explicit in `.github/workflows/ci.yml` and mirrored in a CI-matrix table in the README.
- **Developer tooling (v0.1):** ships `.pre-commit-config.yaml` (ruff for Python lint+format, `clang-format-16` with LLVM style for C++, `cmake-format` for CMakeLists, `check-yaml`, `end-of-file-fixer`), `.editorconfig` (4-space Python, 4-space C++, LF line endings), and `CONTRIBUTING.md` stub. CI runs `pre-commit run --all-files` as a first-check job.

### C++ API shape (v0.1)

Illustrative API sketch for the proposed public surface. File paths in comments are target locations, not current-tree locations. High-level pipelines plus low-level composable primitives all share the same param structs:

```cpp
namespace meshmonk {

    // Canonical types (library/include/meshmonk/types.hpp)
    using FeatureMat  = Eigen::Matrix<float, Eigen::Dynamic, 6>;
    using FacesMat    = Eigen::Matrix<int,   Eigen::Dynamic, 3>;  // signed int32 for numpy/trimesh interop
    using VecDynFloat = Eigen::VectorXf;
    inline constexpr int NUM_FEATURES = 6;

    // Strong rigid transform (library/include/meshmonk/transform.hpp)
    // Semantic contract: maps original floating positions to aligned positions,
    // i.e. result.transform.apply(original_features) == result.aligned_features
    // (up to floating-point rounding).
    // apply(): positions: R*p + t; normals: R*n (rotation only, no translation).
    // inverse(): handles similarity transforms correctly when use_scaling=true
    //   (inverse of [sR|t] is [(1/s)R^T | -(1/s)R^T*t], not just [R^T|-R^T*t]).
    // compose(): matrix multiplication — correct for both SE(3) and similarity.
    struct RigidTransform {
        Eigen::Matrix4f matrix = Eigen::Matrix4f::Identity();   // 4x4 SE(3); similarity if params.use_scaling
        RigidTransform compose(const RigidTransform&) const;    // returns (*this) * rhs
        FeatureMat     apply(const FeatureMat&) const;          // transforms positions (R*p+t); rotates normals (R*n)
        RigidTransform inverse() const;
    };

    // Composable parameter structs (library/include/meshmonk/params.hpp)
    struct CorrespondenceParams { bool symmetric=true; int num_neighbours=3;
                                   float flag_threshold=0.9f; bool equalize_push_pull=false; };
    struct InlierParams { float kappa=12.0f; bool use_orientation=true; };
    struct ViscoElasticParams {
        float sigma = 3.0f;
        // Viscous + elastic annealing pairs interpolate across num_iterations.
        // Defaults mirror nonrigid MATLAB demo (start == 200, end == 1).
        // IMPORTANT: for pyramid registration, the MATLAB demo ties start to
        // num_iterations (= 90). The Python wrapper auto-populates
        // viscous/elastic start = num_iterations when the user does not set
        // them explicitly, matching the MATLAB convention. The C++ struct
        // defaults (200) are fallbacks for the nonrigid case only.
        int num_viscous_iterations_start = 200;
        int num_viscous_iterations_end   = 1;
        int num_elastic_iterations_start = 200;
        int num_elastic_iterations_end   = 1;
    };
    // Percentage of vertices to remove (0.0 = no decimation, 90.0 = remove 90%).
    // Internal code divides by 100 before use. Values match MATLAB demo convention.
    struct DownsampleSchedule { float float_start=50.f, target_start=70.f,
                                      float_end=0.f,   target_end=0.f; };

    struct RigidParams    { CorrespondenceParams correspondences; InlierParams inliers;
                            int num_iterations=80; bool use_scaling=false; };
    struct NonrigidParams { CorrespondenceParams correspondences; InlierParams inliers;
                            ViscoElasticParams transform; int num_iterations=200; };
    struct PyramidParams  { CorrespondenceParams correspondences; InlierParams inliers;
                            ViscoElasticParams transform; DownsampleSchedule downsample;
                            int num_iterations=90; int num_pyramid_layers=3; };

    // Result structs (library/include/meshmonk/result.hpp)
    // v0.1 ships minimal shapes; diagnostic fields (fitness, inlier_rmse,
    // num_inliers, converged) are deferred to v0.2 once detection sites exist.
    //
    // Implementation note: displacement_field in NonrigidResult/PyramidResult is
    // computed at the wrapper level as (aligned_positions - original_positions),
    // NOT extracted from the internal classes (NonrigidRegistration mutates
    // features in-place and destroys the ViscoElasticTransformer's displacement
    // field when the function returns). The wrapper saves a copy of input
    // positions before calling update(), then diffs — cost: one (N,3) float32
    // copy (~650KB for 54K vertices). Similarly, final_inlier_weights requires
    // running a final correspondence + inlier-detection pass after registration
    // completes, since the per-iteration weights are not accumulated.
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
    [[nodiscard]] tl::expected<RigidResult, RegistrationError>
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
    //
    // Lifetime requirement: all Eigen::Ref<const ...> parameters must remain
    // valid for the duration of the function call. Internal KD-tree construction
    // holds references to input data. Do not pass temporaries that would be
    // destroyed before the call returns. (nanobind GC-pins numpy arrays
    // automatically, so Python callers are safe.)
}
```

**Param struct field-order stability rule:** field order inside each param struct is STABLE across minor versions (v0.x → v0.x+1). Reordering fields requires an ADR update. The Python shim uses named `nb::arg()` bindings (not positional construction) so Python callers are immune to field reordering, but C++ callers using designated initializers would be. This protects against silent field-swap hazards (e.g., `num_iterations` and `num_neighbours` both being `int` and adjacent).

**Error detection sites (v0.1):** `RegistrationError` is a real enum, not aspirational. Current code writes `std::cerr` on failure and continues with garbage; the rewrite converts these into expected-style returns at four call sites:

| Error | Detected at | Criterion |
|---|---|---|
| `DegenerateInput` | boundary pre-check | empty inputs, dim/shape mismatch, `FeatureMat` rows ≠ `floating_flags.size()`, zero-range bbox, or all-zero inlier flags |
| `InsufficientInliers` | after `InlierDetector` (boundary-converted) | <4 non-zero inlier weights (rigid needs ≥3 non-collinear points for the cross-variance matrix to yield a unique rotation quaternion) |
| `DecompositionFailed` | internal (boundary-converted) | `RigidTransformer` eigendecomposition (Horn method via `SelfAdjointEigenSolver` on the 4×4 Q matrix) or `EigenVectorDecomposer` fails to converge on valid-shaped input |
| `NonConvergence` | wrapper loop (v0.2+ only) | convergence criterion (TBD) not met within `num_iterations`; **not raised in v0.1** since convergence tracking is deferred |

**Split between boundary and internal detection:** `DegenerateInput` is raised by an explicit pre-check at the top of each public API function (`rigid_registration()`, `nonrigid_registration()`, `pyramid_registration()`) — so shape / rank / dim errors never reach the algorithmic core. Internal algorithmic classes (`ViscoElasticTransformer`, `InlierDetector`, `RigidTransformer`) are **NOT** rewritten to return an expected type; they retain their existing control flow and signal failure via the log sink plus an internal failure-status bridge that the public wrapper converts before returning. This preserves the "load-bearing code untouched" constraint while still giving callers typed errors.

**Portability note:** v0.1 uses `tl::expected` unconditionally (header-only, single file, API-compatible for our use) to sidestep MSVC/libstdc++/libc++ `<expected>` gaps. `std::expected` remains the medium-term target once the compiler matrix is proven.

All internal `std::cerr` writes are routed through a logger sink (configurable via `meshmonk.set_log_level()` in Python) so consumers can suppress stderr noise. The log sink is thread-safe (atomic level read) but its effect is global to the process.

Load-bearing algorithmic code (`ViscoElasticTransformer`, `InlierDetector`, `RigidTransformer`, `NonrigidRegistration` annealing, `PyramidNonrigidRegistration` scheduling) remains untouched. The internal `namespace registration` and `global.hpp` are preserved as-is; the new `namespace meshmonk` public API coexists with it. The redesign is entirely at the API surface above them. The `extern "C"` block and all `_mex` function variants in `meshmonk.hpp` are deleted as part of the v0.1 API redesign (consistent with ADR-001 D3, drop MATLAB).

**"Confine OpenMesh to I/O boundary" — scoped clarification:** this goal applies to *new* public-API-level code. `ViscoElasticTransformer` retains its internal `TriMesh` member for per-iteration normal recomputation (this is load-bearing algorithmic usage, not I/O). `Downsampler::update()` rebuilds a TriMesh per call for quadric decimation — also kept. The promise means: no *new* code reconstructs TriMesh objects; the public API wrapper operates purely on Eigen matrices. Writing a pure Eigen-based `compute_normals` (~30 LOC, face-area-weighted averaging) to replace the OpenMesh path in `ViscoElasticTransformer` is a reasonable v0.1 stretch goal but not required.

**Low-risk internal cleanups (behavioral no-ops):**
- Replace raw `new`/`delete` of `CorrespondenceFilter` in `RigidRegistration::update()` and `NonrigidRegistration::update()` with `std::unique_ptr` (two lines each) — exception-safety fix, no algorithmic change.
- Remove gratuitous OpenMesh `#include` from `ScaleShifter.hpp` — the class uses only Eigen matrices and `NeighbourFinder`; the include is dead weight.

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
# Note: compute_correspondences operates in 6D (positions + normals) —
# inputs must be (N,6) feature matrices, not (N,3) position matrices.
feat_float = meshmonk.features_from_vertices(V_float, F_float)  # (N,6) helper
feat_target = meshmonk.features_from_vertices(V_target, F_target)
corrs = meshmonk.compute_correspondences(feat_float, feat_target, flags_float, flags_target)
weights = meshmonk.compute_inlier_weights(feat_float, corrs, kappa=12.0)
transform = meshmonk.compute_rigid_transform(feat_float, corrs, weights)
normals = meshmonk.compute_normals(V_float, F_float)       # (N,3) vertex normals
V_down, F_down = meshmonk.downsample_mesh(V_float, F_float, ratio=50.0)
V_shifted = meshmonk.scale_shift_mesh(V_float, V_target)   # center + scale to target bbox

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

**nanobind boundary-copy caveats:** the zero-copy path requires `float32` for `FeatureMat` and `int32` for `FacesMat`, plus C-contiguous column-major memory. In practice, numpy arrays are row-major by default, so nanobind performs **one implicit copy per call at the Python→C++ boundary** to match Eigen's column-major storage. For a (54K, 6) float32 matrix this is ~1.3 MB — trivial. Default `trimesh` arrays use `int64` faces — the wrapper does a one-time cast at the Python boundary so end users aren't penalized. If profiling flags the boundary copy as hot in v0.2, internal storage switches to row-major — a larger refactor deferred.

**Normals handling:** precedence is (a) explicit `normals=V_normals` kwarg if passed; else (b) `mesh.vertex_normals` if the mesh object exposes the attribute and its L2 norm per row is non-zero; else (c) recomputed via `meshmonk.compute_normals(V, F)`. Pass `compute_normals=True` to force recomputation even when `.vertex_normals` is available. For the raw-array path, if `floating_features` columns 3–5 are all-zero (as current MATLAB demos pass — `single(zeros(size(floatingPoints)))`), the wrapper auto-recomputes normals and emits a one-line warning via the log sink, matching current C++ behavior while surfacing the smell for the migration. `meshmonk.features_from_vertices(V, F)` is the recommended helper to avoid the warning.

Errors bubble as Python exceptions (nanobind translates unexpected return values via a registered translator). Type stubs (`.pyi`) shipped with the wheel from v0.2 onward.

**GIL release (v0.1 deliverable):** all registration functions and composable primitives release the Python GIL via `nb::call_guard<nb::gil_scoped_release>()` before entering C++. Without this, a `meshmonk.rigid_register()` call (30–180 s) blocks all Python threads for the entire duration. nanobind GC-pins numpy arrays for the duration of the call, so data lifetime is safe under GIL release.

### Python package dependencies

Hard runtime deps (in `[project] dependencies`): `numpy >= 1.24`, `typer >= 0.9`. Build-time requires (`[build-system] requires`): `scikit-build-core >= 0.8`, `nanobind >= 2.0`, `cmake >= 3.26`.

Optional I/O extras (in `[project.optional-dependencies]`):

- `meshmonk[io]` — pulls `trimesh >= 4.0` (chosen over `meshio` for richer mesh-object model, built-in normals computation, and wider format coverage). Used by the wrapper's path-input handling and by `meshmonk demo`.
- `meshmonk[dev]` — pulls `pytest`, `pre-commit`, `ruff` for contributors.

The CLI (`meshmonk demo`, `meshmonk rigid`, etc.) depends on the `io` extra since it loads OBJ files. The `meshmonk` entry point emits a clear `pip install 'meshmonk[io]'` hint if `trimesh` is missing.

### Canonical parameter defaults

Use the MATLAB demo-script values as the Python API defaults. These are grounded in the current repo's `demo/test_rigid_registration.m`, `demo/test_nonrigid_registration.m`, and `demo/test_pyramid_registration.m`, and they differ materially from the current C++ header defaults. Key differences:

| Param | Current C++ | MATLAB (adopted) |
|---|---|---|
| `rigid.num_iterations` | 20 | **80** |
| `inlier_kappa` | 4.0 | **12.0** |
| `correspondences_num_neighbours` | 5 | **3** |
| `correspondences_flag_threshold` (rigid/nonrigid) | 0.99 | **0.9** |
| `nonrigid.num_iterations` | 60 | **200** |
| `correspondences_flag_threshold` (pyramid) | 0.99 | **0.999** |
| `pyramid.num_iterations` | 60 | **90** |

**Param passthrough rule:** the new public API wrapper functions MUST pass *all* params from the param structs through to the internal classes explicitly. Internal class defaults (e.g., `downsampleFloatStart = 90.0f` in `PyramidNonrigidRegistration`, `numIterations = 10` in `NonrigidRegistration`) differ from the adopted MATLAB-derived defaults. The wrapper must always call `set_parameters()` on internal classes before `update()`, never relying on internal member-initializer defaults.

### Build system migration (v0.1 deliverable)

The current `pyproject.toml` uses `setuptools>=77` (v0.0 scaffolding). v0.1 migrates to `scikit-build-core` as the PEP 517 backend — this gates nanobind binding compilation and should be done early.

The current `CMakeLists.txt` uses `FetchContent_Declare` to download Eigen from GitLab at build time, contradicting the overview design's "no FetchContent at build time — `pip install` works offline" rule. v0.1 vendors Eigen 3.4.0 headers under `vendor/eigen-3.4.0/` and updates CMakeLists.txt to use `add_subdirectory()` on the vendored copy. (`vendor/` currently contains only `OpenMesh-11.0.0` and `nanoflann.hpp`.)

### Harness strategy

Tiered harness, ordered fast → slow:

- **Tier 1 — Analytical** (deterministic, math-grounded, <5 s total): synthetic rigid recovery (apply known SE(3) to a mesh's full (N,6) feature matrix — positions AND normals — via `RigidTransform::apply()`, run `rigid_registration` from identity, verify both position and normal recovery), self-consistency, round-trip, correspondence sanity, inlier-weight sanity. The synthetic mesh MUST have non-degenerate normals (not all-zero) since the correspondence search operates in 6D feature space; all-zero normals would degenerate the search to effectively 3D. This is where end-to-end ICP gets covered.
- **Tier 2 — Primitive fixture** (deterministic, <1 s): 8-vertex cube + known correspondences + inlier weights. Exercises `compute_rigid_transform` as a least-squares primitive, **not** end-to-end ICP. Ported from `origin/add-compute-rigid-transform-cli:cli/test_data/rigid_transform/`.
- **Tier 3 — Human-approved visual golden** (captured once per registration type in v0.1): owner inspects output in Blender or MeshLab, approves, we freeze as `.npz`. Future runs compare with Hausdorff/RMSE via `tests/utils/mesh_compare.py` (~20 lines of numpy).
- **Tier 3.5 — Legacy scientific-equivalence gate**: reference outputs on `Template.obj ↔ demoFace.obj` anchor the "did the rewrite preserve the current library's numerical behavior?" baseline. New-code output must match within calibrated tolerance. Trust comes from the *library* (research-validated per ADR D1), not from whatever plumbing wraps it. **Capture source, in priority order:** (a) **checked-in outputs on `origin/cli:demo/`** — `rigid_output.obj`, `pyramid_output.obj`, `rigid_transform.txt` from a prior `meshmonk_cli` run are already committed; v0.0 ingests these directly via `obj_to_npz.py`; (b) reproducing those outputs from a fresh `origin/cli` build (`cmake --build build --parallel` after `cmake -S .worktrees/cli -B /tmp/meshmonk-cli-build`) on the owner's toolchain as a separate **reproducibility check** — divergence is information, not a blocker; (c) the MATLAB demos (`test_rigid_registration.m`, `test_pyramid_registration.m`) as a cross-check via MEX against the same research-validated library. Not circular with ADR D7's rejection of untrusted pybind11 output — the baseline is the raw library, not the abandoned binding layer. **Trust-promotion rule:** Tier 3.5 starts as an **advisory** gate (diff reported, not CI-failing) because the committed goldens and the CLI that produced them landed in a single unverified commit (`origin/cli@3b72b72`, 2025-05-29). Promote to a hard CI gate only once at least one of (b) or (c) has confirmed the goldens independently.
- **Tier 4 — Soft anchor**: `origin/cli:demo/rigid_transform.txt` (4×4 rigid transform from a `meshmonk_cli rigid_reg` run, committed on the `cli` branch) carried into `data/` during v0.0 — loose-tolerance sanity
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
