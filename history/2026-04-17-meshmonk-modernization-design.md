# Design: MeshMonk Modernization — Python-First Rebuild

**Date:** 2026-04-17
**Status:** Draft (strategic direction settled; implementation details still provisional)
**Decisions:** [ADR-001](../docs/decisions/ADR-001-meshmonk-modernization.md)
**Origin:** Brainstorming session 2026-04-17; prompted by revisiting the May 2025 "Phase 5" spec and finding it too narrow
**Supersedes:** May 2025 `modernization-and-python-bindings` spec (remote-branch document; not present in this checkout)
**Phase docs:** v0.0, v0.1, v0.2, v0.3 (linked below)

---

## Problem

**Reference-frame note:** this document intentionally refers to three different things, and they are easy to blur together unless called out explicitly:

- **Current repo checkout** — what exists in this working tree today (`matlab/`, root `example.cpp`, `demo/`, `tutorial/`, current `meshmonk.hpp`)
- **Remote branch state** — work visible under `origin/*` but not in this checkout (for example `origin/cli`, `origin/add-compute-rigid-transform-cli`, `origin/feat-python-bindings`)
- **Proposed post-overhaul layout** — files and directories that do not exist yet, but are the intended target (`library/`, `python/`, `tests/`, `pyproject.toml`)

Unless stated otherwise, references below are proposals, not claims that those paths already exist locally.

MeshMonk is a C++ 3D mesh registration library (rigid ICP, nonrigid deformable registration, pyramid multi-resolution). Algorithmic core represents ~10 years of research-validated tuning by a university lab. A 2025 modernization effort (Phases 1–4) reorganized the repo and added first-pass pybind11 bindings written by a then-young AI agent.

The current state across the current repo plus active remote branches has real structural problems that go beyond what the May 2025 "Phase 5" plan (demo scripts + README polish) would fix:

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

This section describes the **proposed post-overhaul architecture**, not the current checkout.

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
│    rigid_registration(...)      → expected-style Result/Error│
│    nonrigid_registration(...)   → expected-style Result/Error│
│    pyramid_registration(...)    → expected-style Result/Error│
│    compute_correspondences(...) → Correspondences            │
│    compute_inlier_weights(...)  → VecDynFloat                │
│    compute_rigid_transform(...) → expected-style Result/Error│
│    ...                                                        │
│  Deps: Eigen3, OpenMesh (→ meshoptimizer v0.4+), nanoflann   │
└─────────────────────────────────────────────────────────────┘
```

### Repository layout (proposed after v0.1)

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

**Planned deletions from the current repo / remote branches:** `matlab/`, all `_mex` C++ functions, `test_meshmonk_mexing*`, remote-branch-only `cli/`, old `test_bindings.py` on the bindings branch, and current/branch-local example scratch files such as root `example.cpp`. Branches `origin/feat-python-bindings`, `origin/add-compute-rigid-transform-cli`, `origin/feature/*`, and `origin/feat/meshmonk-cli*` are retired once their useful fixtures are extracted.

### Staged roadmap

**v0.0 — Reality capture (~1–2 days):** establish trust anchors before refactoring. Create the `meshmonk-modernization` branch from `origin/cli` and check it out in the dedicated worktree `.worktrees/meshmonk-modernization/` (see Branch strategy). All subsequent work — v0.0 through v0.1 — happens on this single branch in that worktree. Cherry-pick the 8-vertex cube fixture from `origin/add-compute-rigid-transform-cli:cli/test_data/rigid_transform/` to `tests/fixtures/cube_rigid/` before that branch is retired. Move current demo assets (`demo/Template.obj`, `demo/demoFace.obj`, `demo/{Demo,Template}FaceLandmarks.csv`) to `data/`, plus `origin/cli:demo/rigid_transform.txt` (verified present on the `cli` branch, not in current checkout). **Capture legacy baseline:** the `cli` branch already has committed outputs — `origin/cli:demo/rigid_output.obj`, `origin/cli:demo/pyramid_output.obj`, `origin/cli:demo/rigid_transform.txt` — from a prior run of `meshmonk_cli` on `Template.obj ↔ demoFace.obj`. These become the v0.0 Tier 3.5 reference, converted to `.npz` via `tests/utils/obj_to_npz.py` and stored under `tests/golden/legacy_baseline/`. Reproducing them from a fresh `origin/cli` build on the owner's Mac is a **reproducibility check** (see checklist), not a blocker for v0.0 — `cmake` is not currently installed locally. Write the mesh-compare helper. Set up pytest scaffolding. Commit this design doc + ADR. No library changes.

See [v0.0 phase doc](./2026-04-17-meshmonk-modernization-v0.0-design.md) for details.

**v0.1 — Foundation (~1–2 focused weeks):** C++20 overhaul + nanobind bindings + Python package + typer CLI + scikit-build-core + GitHub Actions CI. `uv pip install .` (non-editable) works on the owner's Mac; editable installs tracked via [astral-sh/uv#14383](https://github.com/astral-sh/uv/issues/14383) — fallback is `uv pip install nanobind scikit-build-core cmake ninja && pip install -e . --no-build-isolation` (build-isolation off means build deps must be in the active env). All required harness tiers passing, including **owner runs one Blender/MeshLab visual sign-off session for rigid / nonrigid / pyramid outputs** to freeze Tier 3 goldens. No `meshmonk::legacy::` hedge — the current public API has no external consumers to preserve (MATLAB/C++-CLI/pybind11 all being dropped).

- **CI matrix (v0.1):** Ubuntu 22.04 + Ubuntu 24.04; macOS 14 (arm64) + macOS 13 (x86_64); Python 3.10, 3.11, 3.12, 3.13; CMake 3.26+ installed via `lukka/get-cmake@latest` at a pinned version; gcc 11+, clang 15+, AppleClang 15+. Matrix is explicit in `.github/workflows/ci.yml` and mirrored in a CI-matrix table in the README.
- **Developer tooling (v0.1):** ships `.pre-commit-config.yaml` (ruff for Python lint+format, `clang-format-16` with LLVM style for C++, `cmake-format` for CMakeLists, `check-yaml`, `end-of-file-fixer`), `.editorconfig` (4-space Python, 4-space C++, LF line endings), and `CONTRIBUTING.md` stub. CI runs `pre-commit run --all-files` as a first-check job.

See [v0.1 phase doc](./2026-04-17-meshmonk-modernization-v0.1-design.md) for details.

**v0.2 — Team-ready:** ship type stubs. README rewrite. Migration guide from MATLAB (side-by-side snippets) at `docs/migration-from-matlab.md`. Windows CI if feasible. cibuildwheel dry-run. Add diagnostic fields (`converged`, `fitness`, `inlier_rmse`, `num_inliers`) to result structs once convergence criteria are pinned down. PR `cli` → `master`; rename `master` → `main`.

See [v0.2 phase doc](./2026-04-17-meshmonk-modernization-v0.2-design.md) for details.

**v0.3 — PyPI-ready:** cibuildwheel matrix (manylinux, macOS x86_64 + arm64, Windows). Test PyPI round-trip. MkDocs-material docs site. Publish to PyPI. `pip install meshmonk` works globally.

See [v0.3 phase doc](./2026-04-17-meshmonk-modernization-v0.3-design.md) for details.

**v0.4+ — Future:** OpenMesh → meshoptimizer; reconsider Meson; MCP server; benchmarks.

### Branch strategy

This section refers to **remote branches under `origin/*` plus proposed new branches**. It is not a statement about what already exists locally in this checkout.

```
master (2021, stale)
  ↑ PR at v0.2 milestone, then rename master → main
  │
cli (integration trunk — current)
  ↑ merge at v0.1 milestone (PR: meshmonk-modernization → cli)
  │
meshmonk-modernization  ← single long-lived branch, forked from origin/cli
                          where ALL v0.0 and v0.1 work happens;
                          v0.0 and v0.1 are tags/milestones on this branch,
                          not separate branches
```

**Development environment:** all modernization work is done inside a dedicated git worktree to keep it isolated from the `master` checkout (which keeps design docs, tmp/ notes, and any browsing):

```
/Users/jonat/code/personal/meshmonk/                       ← master checkout (design, notes)
/Users/jonat/code/personal/meshmonk/.worktrees/
    meshmonk-modernization/                                ← active dev worktree
    cli/                                                   ← read-only reference (origin/cli)
    cube-fixture/                                          ← read-only reference
```

Create it once during v0.0:

```
git worktree add -b meshmonk-modernization .worktrees/meshmonk-modernization origin/cli
cd .worktrees/meshmonk-modernization
```

Editor/IDE sessions and `cmake` builds run from that worktree; design-doc edits and brainstorming continue in the `master` checkout. This avoids constantly swapping branches in a single working tree and keeps `.worktrees/cli` and `.worktrees/cube-fixture` as cheap, read-only references for legacy baseline ingestion and cherry-picking. The `.worktrees/cli` tree is "read-only" in the sense that `origin/cli` is never modified or pushed from there — a disposable `build/` (or `_deps/`) directory produced by the v0.0 reproducibility check is expected and fine; wipe it afterward.

**Tradeoff accepted:** `meshmonk-modernization` is a long-running (~1–2 weeks) unreviewed branch because the work is solo. Mitigation: when the `v0.0` tag is cut, self-review the full diff against `origin/cli` before starting v0.1 refactoring — and once a tag is pushed, treat it as immutable (no force-push, no re-tagging).

Branches to retire: all remote branches except `master`, `cli`, and the new `meshmonk-modernization` — explicitly including `feat-python-bindings`, `add-compute-rigid-transform-cli` (after cube-fixture cherry-pick), `feature/matlab-cmake-adaptation`, `feature/phase1-cmake-restructure`, `feat/meshmonk-cli`, `feat/meshmonk-cli-superbuild-attempt`, `cpp_project`, `development`, `mh-fixdocs`.

### Convention: floating vs. target

**Enforced across entire API:** `floating` is the mesh being transformed (conceptually: the template we deform); `target` is the fixed reference (the face scan we're trying to match). Verified consistent across all MATLAB demos.

---

## Key Design Decisions

All decisions captured in [ADR-001](../docs/decisions/ADR-001-meshmonk-modernization.md):

- **D1 (FIRM):** Modernize in C++20 — not Rust
- **D2 (FIRM):** nanobind for Python bindings — not pybind11
- **D3 (FIRM):** Drop MATLAB support — redirect to university fork
- **D4 (FIRM):** Delete C++ CLI — replace with Python typer CLI
- **D5 (FIRM):** CMake + scikit-build-core — not Meson
- **D6 (FIRM):** Params structs + result structs + expected-style typed failure
- **D7 (FLEXIBLE):** Tiered harness, analytical first, no circular baselining
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
- Branch structure: retirement of 9+ obsolete branches; one new long-lived `meshmonk-modernization` branch (v0.0 and v0.1 land as tags/milestones on it) developed in a dedicated worktree
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
- **PyPI trusted-publishing setup (v0.3 prerequisite):** see [v0.3 phase doc](./2026-04-17-meshmonk-modernization-v0.3-design.md).
- **Repository governance**: the repo was transferred from `TheWebMonks/meshmonk` to `jsnyde0/meshmonk` on 2026-04-17 (preserves the 85 stars; existing `git clone` URLs continue to work via GitHub's 301-redirect; anchors the Python rewrite under personal attribution). PyPI namespace `meshmonk` still to be claimed by `jsnyde0`. If TheWebMonks wants to continue MATLAB work, they fork from the transferred repo.
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

## Phase documents

- [v0.0 — Reality Capture & Scaffolding](./2026-04-17-meshmonk-modernization-v0.0-design.md)
- [v0.1 — Foundation](./2026-04-17-meshmonk-modernization-v0.1-design.md)
- [v0.2 — Team-ready](./2026-04-17-meshmonk-modernization-v0.2-design.md)
- [v0.3 — PyPI-ready](./2026-04-17-meshmonk-modernization-v0.3-design.md)
