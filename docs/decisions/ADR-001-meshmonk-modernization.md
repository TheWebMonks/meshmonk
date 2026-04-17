# ADR-001: MeshMonk Modernization — Python-First Rebuild

**Status:** Accepted
**Date:** 2026-04-17
**Design:** [Design doc](../../history/2026-04-17-meshmonk-modernization-design.md)
**Parent:** —
**Related:** Supersedes the May 2025 `modernization-and-python-bindings` spec (remote-branch document; not present in this checkout)

**Scope note:** this ADR captures strategic decisions. File paths and branch names cited in the design doc are expected to be labeled as one of: current checkout, remote branch state under `origin/*`, or proposed future layout.

## Context

MeshMonk is a C++ 3D mesh registration library developed by a university research group (~10 years of algorithmic tuning). Between May and June 2025 a modernization effort (Phases 1–4) restructured the repo, introduced a CMake build system, and produced a first-pass pybind11 Python binding (`meshmonk_python`) — largely written by an older AI coding agent with limited validation.

The 2025 plan's Phase 5 was scoped as "demo/test scripts + README polish." Revisiting that plan in 2026 revealed deeper issues that make Phase 5 the wrong frame:

- The pybind11 bindings are untested and untrusted (author's own comments note confusion about the API shape).
- There is no `pyproject.toml`, no `pip install` story, no CI.
- The C++ public API has 15+ positional-arg telescoping signatures, in-place mutation, and no typed failure mode.
- ~45% of the C++ codebase is MATLAB MEX plumbing that blocks build-system flexibility and adds maintenance burden.
- The C++ CLI was built by the same untrusted AI-agent run and was never validated.
- The repo has essentially no automated test coverage — no pytest, no integration tests, no CI.

The user explicitly wants the "cleanest best architecture to be proud of," is open to radical rewrites (including Rust), and has no sunk-cost attachment to the existing code. This ADR captures the strategic decisions that shape the resulting overhaul.

## Decisions

### D1: Modernize in C++20 — do not rewrite in Rust

**Firmness: FIRM**

Keep the algorithmic core in C++, upgrade to C++20, and redesign the public API surface. Do not port to Rust.

**Rationale:**

The load-bearing algorithmic value of MeshMonk — viscoelastic diffusion, SVD-based rigid fit, Mahalanobis inlier detection, pyramid annealing — represents roughly 30% of the codebase and reflects years of numerical tuning. This is exactly the code most expensive and most risky to rewrite.

The Rust scientific-geometry ecosystem is not mature enough in 2026 to support this work credibly:
- No production-quality halfedge mesh library (what we use OpenMesh for). `plexus` is unmaintained.
- `nalgebra`/`faer` cover basic linear algebra but lack Eigen's ecosystem depth for SVD/sparse solvers.
- No competitor in this space (Open3D, libigl, PyMeshLab, CGAL) is Rust-based.

Realistic cost of a Rust rewrite with numerical re-validation: 12–18 months. The memory-safety argument does not apply to a batch scientific tool processing trusted input.

C++20 gives us the modern idioms we actually need: designated initializers for aggregate structs, `[[nodiscard]]`, concepts, and a clean path to standard-library `expected` once the compiler matrix permits it.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **Modernize in C++20 (chosen)** | • Preserve algorithmic value<br>• Eigen/OpenMesh ecosystem intact<br>• Clean idioms available<br>• Weeks of work, not years | • Still C++ (verbose, build complexity)<br>• Manual memory management risks |
| Rust rewrite (PyO3) | • Memory safety<br>• Modern tooling<br>• Smaller wheels | • No OpenMesh equivalent<br>• Ecosystem immature for geometry<br>• 12–18 month revalidation<br>• No field precedent |
| Stay on C++14, narrow Python scope | • Minimal risk | • Keeps all current API flaws<br>• No path to modern idioms |

**What would invalidate this:**

- A production-quality Rust halfedge mesh library appears and is adopted by 2+ major geometry projects
- The algorithmic code is fully re-derived and re-verified against a new numerical reference, making the C++ tuning non-unique
- Performance profiling reveals a fundamental C++ limitation that Rust solves

### D2: nanobind for Python bindings — not pybind11

**Firmness: FIRM**

Replace the existing pybind11 bindings with nanobind. Treat the existing pybind11 wrapper as disposable.

**Rationale:**

nanobind is the pybind11 author's (Wenzel Jakob) rewritten successor, actively developed and used in production by JAX, Mitsuba 3, DrJit, and — decisively — libigl-python-bindings, the closest analogue to what MeshMonk-Python should be.

Concrete wins for this project:
- 2–10× smaller compiled binary (Eigen-heavy extensions go from ~800 KB to ~120–200 KB)
- 2–5× faster compile times
- 2–3× faster cold import
- First-class `stubgen` for shipping `.pyi` type stubs
- Near-identical Eigen bridge API (`nanobind/eigen/dense.h` ≈ `pybind11/eigen.h`)

Migration cost is a single focused afternoon. The 330-line Jules-written pybind11 file has no sunk cost to protect.

**Python version + ABI strategy:**

- v0.1 targets **CPython 3.10–3.13**. Older Pythons are EoL or close to it; 3.10 gives `match`, built-in generic types (`list[int]`), and `ParamSpec` without `from __future__ import annotations`.
- v0.3 wheels are built with nanobind's **stable-ABI mode** (`NB_STABLE_ABI=1`), collapsing the wheel matrix from ~4 Pythons × 3 OSes × 2 arches = 24 wheels to 6 wheels (one per OS/arch). abi3 has <1% typical overhead; if profiling later shows it's >2%, flip the build flag for v0.4 — no API impact.
- Build-time numpy is pinned to **2.0** for forward-compatibility with numpy 1.x at runtime.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **nanobind (chosen)** | • Smaller wheels<br>• Faster builds<br>• Built-in stubgen<br>• libigl precedent | • C++17 minimum (we need C++20 anyway)<br>• Less Stack Overflow history |
| Continue pybind11 | • Familiar, vast ecosystem | • Larger binaries<br>• Slower builds<br>• External stubgen tools |
| Cython | • Mature | • Poor C++ binding ergonomics |
| SWIG | • Legacy | • Verbose, hard to maintain |
| ctypes/cffi | • Minimal binary size | • No type safety, no NumPy integration |

**What would invalidate this:**

- nanobind's maintenance pauses or Wenzel Jakob steps back
- A blocker emerges in nanobind's Eigen/NumPy interop that doesn't affect pybind11
- Binary size regresses to pybind11 parity

### D3: Drop MATLAB support

**Firmness: FIRM**

Remove all MATLAB MEX code (`matlab/`, `_mex` function variants, and MATLAB-specific demo/tooling). Redirect MATLAB users to the university fork that already exists.

**Rationale:**

MATLAB support is ~45% of the non-algorithmic code in the current codebase. It requires:
- Maintaining parallel `_mex` raw-pointer API surfaces
- Keeping CMake as the build system (Meson has no MEX support — though CMake is independently the right call, see D5)
- MATLAB-specific demo tooling
- Cross-platform MEX testing that has never existed

Users who need MATLAB can use the university fork. This is explicitly sanctioned by the user (repo owner) as the official MATLAB path going forward. A README banner points MATLAB users there.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **Drop MATLAB (chosen)** | • ~250 LOC deleted<br>• Single API surface<br>• Simpler tests, CI, release | • MATLAB users must migrate to fork |
| Keep MATLAB minimally | • Backward compatibility | • Perpetual parallel-API maintenance<br>• Blocks clean redesign |
| Deprecate over multiple releases | • Gentle migration | • Slows down primary Python work |

**What would invalidate this:**

- The university fork becomes unavailable or unmaintained
- A large paying customer specifically needs MATLAB support in this repo
- Python bindings fail to reach MATLAB-user feature parity within 2 minor versions

### D4: Delete the C++ CLI, replace with a Python `typer`-based CLI

**Firmness: FIRM**

Delete the entire `cli/` directory on the CLI branch (`cli.cpp`, `cxxopts` dep, shell test scripts). Build the CLI in Python using `typer`, distributed as a `meshmonk` entry point via `pyproject.toml`.

**Rationale:**

The current C++ CLI was built by an older AI agent, has no numerical validation, and was under active churn (see branch `add-compute-rigid-transform-cli`, which investigation confirmed is 100% plumbing with zero library-level value). A research agent verified the branch adds no unique math — every CLI command just calls a library function that already existed.

A Python CLI is materially simpler:
- ~100 lines of Python vs. ~600 lines of C++
- Uses `meshio`/`trimesh` for I/O — format-agnostic (`.obj`, `.ply`, `.stl`, `.vtk`) for free
- Trivially testable with pytest
- Agent-friendly (typed arguments, auto-help, shell completion)
- Share types and behavior with the Python library layer

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **Python CLI via typer (chosen)** | • Simpler<br>• Format-agnostic I/O<br>• Trivial to test<br>• Shared with library | • Requires Python at runtime |
| Keep C++ CLI, finish Karlijne's verification | • Binary-only deployment | • Pure plumbing, no trust, no tests<br>• Dead-end path |
| Rust CLI via clap | • Standalone binary | • New language in build<br>• No benefit over Python |

**What would invalidate this:**

- A user materializes who needs the CLI without Python runtime (embedded, air-gapped), in numbers justifying a separate binary path
- Python CLI perf becomes a bottleneck

### D5: CMake + scikit-build-core — do not migrate to Meson

**Firmness: FIRM**

Keep CMake as the primary build system. Use `scikit-build-core` as the PEP 517 backend for Python packaging. Don't migrate to Meson.

**Rationale:**

OpenMesh is a CMake project. As long as we keep it in-tree for v0.1 (see D8), CMake must be in the toolchain. Meson would require either `cmake.subproject()` layering (more complex than pure CMake) or porting OpenMesh to Meson (multi-month commitment).

Beyond the OpenMesh constraint, scikit-build-core is independently the correct 2026 answer for CMake-based Python packages:
- PEP 517 native (no `setup.py`)
- First-class `cibuildwheel` support (for v0.3 PyPI)
- Minimal changes to existing `CMakeLists.txt` — just add `pyproject.toml`
- Used by libigl-python, Awkward Array, PyBaMM, and nanobind's own examples

Known rough edge: `uv pip install -e .` (editable) has an [open bug](https://github.com/astral-sh/uv/issues/14383) with scikit-build-core. Workaround: `uv pip install .` (non-editable) for production, `pip install -e . --no-build-isolation` for development.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **CMake + scikit-build-core (chosen)** | • OpenMesh-compatible<br>• PEP 517 modern<br>• cibuildwheel ready<br>• libigl precedent | • uv editable-install rough edge (workaround available) |
| Meson + meson-python | • Cleaner DSL<br>• NumPy/SciPy adopted | • Requires `cmake.subproject()` for OpenMesh<br>• Or port OpenMesh (months)<br>• No net gain |
| setuptools + CMakeExtension | • Legacy familiarity | • Deprecated pattern<br>• Editable installs unreliable |

**What would invalidate this:**

- We drop OpenMesh entirely (see D9 — post-v0.4+) AND measure real pain from CMake complexity
- scikit-build-core stalls as a project
- A new build backend emerges that strictly dominates for CMake projects

### D6: Redesign the public API — params structs + result structs + expected-style typed failure

**Firmness: FIRM**

Replace the current telescoping 15-arg in-place-mutation C++ API with:
- Config **aggregate structs** (`RigidParams`, `NonrigidParams`, `PyramidParams`) using C++20 designated initializers
- Value-returned **result structs** (`RigidResult`, `NonrigidResult`, `PyramidResult`)
- `[[nodiscard]]` expected-style return values for typed failure (`tl::expected` in v0.1; `std::expected` once the compiler matrix permits it)
- `Eigen::Ref<const T>` for zero-copy input parameters
- Strong-typed `RigidTransform` (with `.compose()`, `.apply()`, `.inverse()`)
- Free functions, no class-based builders

Expose both **high-level pipelines** (`rigid_registration`) and **low-level primitives** (`compute_correspondences`, `compute_inlier_weights`, `compute_rigid_transform`) sharing the same param structs.

The Python API layers over this with kwargs-friendly shims and duck-typed mesh objects (`.vertices`, `.faces`, `.flags` attributes).

**Rationale:**

Derived from first principles, not library cargo-culting:
- **Minimize numpy-to-C++ friction:** `Eigen::Ref<const T>` is zero-copy for *read-only* inputs from `Eigen::Map` (what nanobind produces from numpy arrays). The floating-features buffer is mutated in place by the internal algorithm, so the wrapper copies it once per call — still a single copy of ~170 KB for the default Template mesh.
- **Introspectable output:** result structs expose what exists today. v0.1 ships minimal fields. `RigidResult` holds `aligned_features`, `transform`, `iterations_run`. `NonrigidResult` and `PyramidResult` additionally expose `displacement_field` (`(N,3)` per-vertex displacement from original floating) so MATLAB-migrants can recover the transform without diffing outputs. Diagnostic fields (`converged`, `fitness`, `inlier_rmse`, `num_inliers`) land in v0.2 once convergence criteria and inlier thresholds are pinned down. Shipping fields with undefined semantics would be worse than omitting them.
- **Composable:** Advanced users run stages individually; primitives share types
- **Safe by construction:** Designated initializers eliminate positional-arg transposition. Expected-style returns make failure explicit for four concrete detection sites: `DegenerateInput` (empty inputs, dim/shape mismatch, zero-range bbox, all-zero inlier flags), `InsufficientInliers` (<4 non-zero inlier weights — rigid needs ≥3 for SVD well-posedness), `DecompositionFailed` (`RigidTransformer` SVD or `EigenVectorDecomposer` fails to converge on valid-shaped input), `NonConvergence` (v0.2+ only; v0.1 runs fixed iterations). Current `std::cerr` failure writes are routed through a logger sink consumers can suppress. Strong types prevent `Transform * Features` confusion.
- **Boundary vs. internal detection split:** `DegenerateInput` is raised by an explicit pre-check at the top of each public API function so shape / rank / dim errors never reach the algorithmic core. Internal algorithmic classes (`ViscoElasticTransformer`, `InlierDetector`, `RigidTransformer`) are **NOT** rewritten to return an expected type — they retain existing control flow and signal failure via the log sink plus an internal failure-status bridge that the public wrapper converts before returning. This preserves the "load-bearing code untouched" constraint while still giving callers typed errors.
- **Pythonic at the top:** kwargs + dataclass-like results + exceptions (translated from `expected`)
- **Field-order stability:** param struct field order is STABLE across minor versions (v0.x → v0.x+1); reordering requires an ADR update. The nanobind shim uses named `nb::arg()` bindings so Python callers are immune to field reordering. C++ callers using designated initializers would still fail to compile on reorder.

**Portability note:** v0.1 uses `tl::expected` unconditionally (header-only, single file, identical API) to sidestep MSVC/libstdc++/libc++ `<expected>` version gaps. v0.2 migrates to `std::expected` once the compiler matrix is known.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **Params struct + result struct + expected (chosen)** | • Discoverable<br>• Safe<br>• Zero-copy<br>• Composable | • C++20 required |
| Fluent builder (`.with_iters(40).run()`) | • Familiar from Java | • Harder to bind to Python<br>• Each setter needs a binding |
| Keep telescoping 15-arg + in-place | • Minimal change | • Every current flaw retained |
| Heavy OO with `Mesh`/`Registration` classes | • Some "OO cleanliness" | • Extra indirection<br>• Poor numpy interop<br>• API bloat |

**What would invalidate this:**

- C++20 designated initializers break on a required target compiler
- Benchmarks show `Eigen::Ref` forces hidden copies for our common call patterns beyond the expected floating-buffer mutation copy

### D7: Harness strategy — tiered, analytical+constructed first, no circular baselining

**Firmness: FLEXIBLE**

Design verification around:
1. **Tier 1 — Analytical** (synthetic rigid recovery, self-consistency, round-trip, correspondence sanity) — grounded in math, no reference needed. End-to-end rigid ICP is exercised here (apply known SE(3), recover from identity), not in Tier 2.
2. **Tier 2 — Primitive fixture** (port the 8-vertex cube fixture from `origin/add-compute-rigid-transform-cli:cli/test_data/rigid_transform/`) — exercises `compute_rigid_transform` as a least-squares primitive, not end-to-end ICP
3. **Tier 3 — Human-approved visual golden** (user inspects in Blender/MeshLab once, we freeze)
4. **Tier 3.5 — Legacy scientific-equivalence gate** — lock reference outputs on `Template.obj ↔ demoFace.obj` as the v0.1 numerical reference. Capture source in priority order: (a) **checked-in outputs on `origin/cli:demo/`** (`rigid_output.obj`, `pyramid_output.obj`, `rigid_transform.txt`) from a prior `meshmonk_cli` run — v0.0 ingests these directly, no build required; (b) reproducing them from a fresh `origin/cli` build on the owner's toolchain as a separate reproducibility check (divergence logged, not a blocker); (c) MATLAB demos (`test_rigid_registration.m`, `test_pyramid_registration.m`) as cross-check or fallback — same research-validated library via MEX. All Tier 3 and 3.5 goldens stored as `.npz` (compressed ~1.3 MB each for 54K-vertex meshes vs 10+ MB as OBJs). NOT circular with rejecting pybind11 output: the baseline is the raw library, not the abandoned binding layer.
5. **Tier 4 — Soft anchor** — `origin/cli:demo/rigid_transform.txt` (committed 4×4 rigid transform from `meshmonk_cli rigid_reg`) carried into `data/` during v0.0 — loose tolerance
6. **Tier 5 — E2E smoke** (CLI commands exit 0, output exists, is valid OBJ)

Explicitly reject freezing **untrusted pybind11 output** as "ground truth" (circular). The Tier 3.5 legacy baseline is different — it captures the raw C++ library behavior that the research group has been using for years.

Visualization is outsourced to Blender / MeshLab. A ~20-line numpy Hausdorff/RMSE helper covers automated golden comparison — no custom viewer.

**Rationale:**

The repo has essentially no existing test coverage. The only mechanically-verified artifact is the 8-vertex cube fixture on `add-compute-rigid-transform-cli`. MATLAB demos represent years of researcher use but have no automated comparison mechanism.

The right move is to start with what math guarantees (Tier 1), add the one existing trustworthy fixture (Tier 2), get one human sign-off per registration type (Tier 3), and build up from there. This avoids the anti-pattern of interpreting code correctness as equivalent to numerical correctness.

Flexibility: the tier selection and tolerance thresholds will evolve in v0.1 as we measure run-to-run numerical noise and identify gaps.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **Tiered harness with analytical + human approval (chosen)** | • No circular trust<br>• Fast CI for Tier 1/2<br>• Human review catches subtle issues | • Requires 1 human session per golden<br>• Tolerance calibration takes measurement |
| Freeze current pybind11 output as golden | • Fast to set up | • Circular (trust untrusted code)<br>• Locks in any existing bugs |
| Only synthetic tests, no real-data goldens | • Fully automated | • Misses algorithm-to-implementation gaps<br>• No real-face coverage |
| Manual testing only | • Simple | • No regression safety net<br>• Blocks CI |

**What would invalidate this:**

- A trustworthy external baseline emerges (published reference implementation with golden outputs)
- Human visual inspection reveals we can't reliably distinguish "right" from "slightly wrong" for nonrigid — would need image-based metrics instead
- CI time for Tier 3 exceeds budget — move Tier 3 to nightly, keep Tier 1/2 per-push

### D8: Keep OpenMesh for v0.1; migrate to `meshoptimizer` in v0.4+

**Firmness: FLEXIBLE**

In v0.1, preserve OpenMesh as the halfedge / quadric-decimation dependency, but confine it to an I/O boundary (no `TriMesh` reconstruction inside `update()` calls).

In v0.4+, evaluate migrating away from OpenMesh. `meshoptimizer` (MIT, Arseny Kapoulkine) covers **only** the quadric-decimation half of what we use OpenMesh for — it is a mesh-data-layout library, not a half-edge topology library. A replacement path requires:
- `meshopt_simplify` for quadric decimation (near-equivalent to `OpenMesh::Decimater::ModQuadricT`)
- **Plus a small half-edge-lite module** (~150–250 LOC) providing boundary-edge detection and vertex-neighbour queries for the `Downsampler`'s boundary-vertex locking. Reference implementations: libigl's `igl::is_border_vertex`, trimesh's `edges_unique`.

Because this is NOT a like-for-like dependency swap, the v0.4+ migration deserves its own design pass, not a one-line roadmap bullet.

**Rationale:**

Replacing OpenMesh in v0.1 would add 1–2 weeks of scope with numerical re-validation risk for the quadric decimation in `Downsampler`. Not worth it for the initial overhaul.

Deferring the migration also keeps the door open for reconsidering Meson (see D5 "what would invalidate").

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **Keep in v0.1, defer migration (chosen)** | • Limits scope<br>• Preserves validated decimation math | • Extra dep, extra build time, bigger binary |
| Drop OpenMesh in v0.1 | • Simpler deps immediately | • Extra weeks of work, re-validation risk |
| Keep OpenMesh forever | • No migration ever | • Misses future simplification opportunity |

**What would invalidate this:**

- OpenMesh development stops or breaks on a needed platform
- Wheel size becomes a real user complaint
- A compelling reason to move to Meson materializes (needs OpenMesh drop first)

### D9: Transfer repo to `jsnyde0/meshmonk`; claim PyPI `meshmonk`

**Firmness: FIRM**

Transfer `TheWebMonks/meshmonk` to `jsnyde0/meshmonk` using GitHub's native repo-transfer feature before v0.1 merges. Claim PyPI namespace `meshmonk` under `jsnyde0`.

**Rationale:**

The repo owner is also the owner of `TheWebMonks`, so this is a unilateral call, not a negotiation. GitHub's repo transfer preserves stars (currently 85), forks, watchers, issues, PRs, releases, and automatically redirects the old URL — existing `git clone` URLs and inbound links keep working. The Python-first rewrite is personal, unpaid work and deserves personal attribution on the owner's GitHub profile. Tying the rewrite to `jsnyde0/meshmonk` also anchors the PyPI namespace naming under the same owner.

If `TheWebMonks` wants to continue MATLAB-flavored work, they fork from the transferred repo — the MATLAB fork then becomes their fork of `jsnyde0/meshmonk`, which is the right topology given the Python-first primary going forward.

**Licensing continuity:** Apache-2.0 continues across the full tree — Apache-2.0 is irrevocable on the 2017 contribution set and this ADR adds nothing to prohibit. A `NOTICE` file at repo root attributes original copyright to TheWebMonks (2017) and new copyright to jsnyde0 (2026). `pyproject.toml` `[project]` lists `license = {text = "Apache-2.0"}`, `authors = [{name = "jsnyde0"}]`. New source files carry `Copyright 2026 jsnyde0` headers; untouched legacy files keep their original `Copyright 2017 WebMonks` headers.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **Transfer to `jsnyde0` (chosen)** | • Preserves 85 stars<br>• Redirects old URLs<br>• Personal attribution<br>• Clean PyPI ownership | • Removes repo from TheWebMonks org page |
| Fork to `jsnyde0` without transfer | • TheWebMonks keeps star count | • Lose star history on rewrite<br>• Two diverging upstreams<br>• MATLAB users confused about canonical location |
| Stay on `TheWebMonks` | • No change | • No personal attribution<br>• Mixes org-scale branding with personal weekend project |

**What would invalidate this:**

- A concrete stakeholder (funder, institution, co-author) requires the repo live on `TheWebMonks` for branding / grant-reporting reasons
- PyPI `meshmonk` is already taken (verify before v0.3; fall back to `meshmonk-py` or similar)

## Related

- [Design doc](../../history/2026-04-17-meshmonk-modernization-design.md) — what we're building and how, implementation details
- Obsolete prior spec: May 2025 `modernization-and-python-bindings` document on a remote modernization branch, not in this checkout
- Research reports dispatched 2026-04-17 (in-session transcripts): pybind11 vs nanobind, Python build backends, mesh library patterns, C++ architecture review, Rust vs C++, modern C++ idioms
