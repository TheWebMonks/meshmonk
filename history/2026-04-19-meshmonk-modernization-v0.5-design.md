# Design: MeshMonk v0.5 — Next Generation

**Date:** 2026-04-19
**Status:** Exploratory (no detailed design — items need individual design passes)
**Parent:** [Overview](./2026-04-17-meshmonk-modernization-design.md)
**Decisions:** [ADR-001](../docs/decisions/ADR-001-meshmonk-modernization.md)

---

## Scope

Architecture-level changes deferred from the v0.1–v0.4 roadmap. These items are tightly coupled and mostly sequential. This doc captures WHAT is planned and WHY, not HOW — each item needs its own design pass before implementation begins.

The parent design doc summarizes these as a single line: "v0.4+ — Future: OpenMesh → meshoptimizer; reconsider Meson; MCP server; benchmarks." This doc expands on what that line actually means and what still needs to be figured out.

---

## Items

### 1. OpenMesh → meshoptimizer

**Status:** Needs design pass (ADR-001 D8 explicitly requires this)
**ADR reference:** D8 (FLEXIBLE) — "Because this is NOT a like-for-like dependency swap, the v0.4+ migration deserves its own design pass, not a one-line roadmap bullet."

**What we know:**

OpenMesh is currently used in two roles inside `Downsampler`:
1. Quadric decimation — `OpenMesh::Decimater::ModQuadricT` reduces vertex count at each pyramid level
2. Topology queries — boundary-edge detection and vertex-neighbour iteration for boundary-vertex locking (so boundary vertices are not collapsed during decimation)

`meshoptimizer` (MIT, Arseny Kapoulkine) replaces only role 1 via `meshopt_simplify`. It is a mesh-data-layout library, not a halfedge topology library. Role 2 has no equivalent in meshoptimizer.

A replacement path therefore requires two components:
- `meshopt_simplify` for quadric decimation
- A half-edge-lite module (~150–250 LOC) providing boundary-edge detection and vertex-neighbour queries. Reference implementations exist: libigl's `igl::is_border_vertex`, trimesh's `edges_unique`. This is not novel work, but it must be written, tested, and numerically validated — decimation results are sensitive to which vertices are locked.

**What we don't know:**

- Whether `meshopt_simplify`'s decimation quality is numerically equivalent to OpenMesh's `ModQuadricT` on MeshMonk's typical inputs (50K–100K-vertex face meshes). The algorithms are similar in principle but differ in implementation details. The Tier 3.5 legacy-baseline gate from v0.1 becomes the acceptance criterion here, but tolerances may need to be relaxed.
- How much of the `Downsampler` internals depend on OpenMesh's `TriMesh` representation beyond decimation and topology queries. A full audit is needed.
- Whether the half-edge-lite module should be a standalone header, part of `library/include/meshmonk/`, or an inline detail file. This affects the public API surface and future maintenance burden.

**What this enables:**

- Removes the largest binary dependency. OpenMesh 11.0.0 dominates the sdist (~20–25 MB of the ~30 MB total). Dropping it would shrink wheels materially and reduce build times.
- Unlocks Meson migration (see item 2) — D5 explicitly conditions Meson reconsideration on OpenMesh being dropped.

**Open questions:**

1. Do we write the half-edge-lite module from scratch, or thin-wrap a smaller existing library (e.g., halfedge header from libigl, which is LGPL — license incompatibility with Apache-2.0 requires investigation)?
2. What is the acceptable tolerance regression in Tier 3.5 goldens after the decimation swap? Needs empirical measurement.
3. Does the migration affect `PyramidResult.displacement_field` semantics? Pyramid-level mesh counts change if decimation quality differs.

---

### 2. Meson migration

**Status:** Contingent on item 1
**ADR reference:** D5 (FIRM) — "What would invalidate: We drop OpenMesh entirely (see D9 — post-v0.4+) AND measure real pain from CMake complexity."

**Condition for this item to proceed:**

Item 1 must be complete. CMake stays as long as OpenMesh is in-tree. Meson's `cmake.subproject()` layering adds complexity that erases any Meson ergonomics gain. This is not negotiable — it is stated explicitly in D5.

**If OpenMesh is dropped, does Meson still make sense?**

Maybe. The argument for Meson becomes: cleaner DSL, faster configuration, adopted by NumPy and SciPy. The argument against: scikit-build-core works well with CMake already, the remaining dependencies (Eigen, nanoflann, meshoptimizer, nanobind) all have first-class CMake support, and there is no measured pain from the current CMake setup after OpenMesh is gone.

D5 lists two conditions that must both be true before Meson is reconsidered: OpenMesh dropped AND real measured pain from CMake complexity. Dropping OpenMesh satisfies only the first. The second must be evaluated at that time — if `CMakeLists.txt` is clean and builds are fast, there is no reason to migrate.

**What we don't know:**

- Whether `meson-python` (the PEP 517 backend for Meson) has matured enough by v0.5 to match scikit-build-core's `cibuildwheel` integration.
- Whether the remaining C++ codebase, once OpenMesh is gone, presents any CMake pain worth solving.

**If Meson migration proceeds, what changes:**

- `CMakeLists.txt` → `meson.build` at root and in `library/`, `python/`
- `pyproject.toml` backend changes from `scikit-build-core` to `meson-python`
- CI matrix job names change; build invocations change
- No API changes — purely build-system

**This item may never happen.** If item 1 is complete and CMake is working fine, Meson migration would be pure churn with no user-visible benefit. The decision gets made at item 1 completion.

---

### 3. MCP server

**Status:** Undefined
**ADR reference:** Non-goal in parent design doc (line 44) — "MCP server wrapping the library for agents — deferred to v0.4+."

**What an MCP server would provide:**

An agent-facing interface wrapping meshmonk operations as MCP tools. Rather than `import meshmonk` from Python, an agent calls `meshmonk/rigid_registration` via MCP protocol and gets a structured result. This is useful when the agent environment does not have a Python runtime with meshmonk installed, or when the agent orchestration layer is language-agnostic.

**What is completely undefined:**

- Which operations to expose as MCP tools. Candidates: `rigid_registration`, `nonrigid_registration`, `pyramid_registration`, individual primitives. Exposing all primitives may be too granular; exposing only pipelines may be too coarse.
- How mesh data is transferred over MCP. MCP is a text/JSON protocol; mesh data is binary. Options: base64-encoded OBJ/PLY, file-path references (server and client share a filesystem), or a streaming binary extension. Each has major tradeoffs.
- MCP protocol version. The spec was at 2025-03-26 as of this writing; it may have moved.
- Whether the MCP server is a separate package (`meshmonk-mcp`) or bundled with the main package behind an optional dependency group.
- Whether the server is implemented in Python (wrapping the Python package) or as a separate process.
- What authentication and sandboxing look like — mesh registration on untrusted input is not a security boundary, but resource limits (mesh size, timeout) are still needed.

**Dependency chain:** MCP server can be implemented independently of items 1 and 2. It wraps whatever version of the Python package is current. However, it makes more sense to do this after the Python API stabilizes (post-v0.2 or v0.3) so the MCP tool signatures don't need to change.

---

### 4. Benchmarks

**Status:** Undefined
**ADR reference:** "benchmarks" appears in the parent design doc's v0.4+ line only. No further specification exists.

**What needs benchmarking:**

At minimum, the performance regression suite needs to cover:
- Registration wall-clock time (rigid, nonrigid, pyramid) on the `Template.obj ↔ demoFace.obj` pair — the Tier 3.5 reference meshes
- Peak memory usage during registration
- Pyramid convergence behavior (iterations-to-threshold) as a proxy for algorithm stability

Potentially also:
- Import time (`import meshmonk` cold start)
- Wheel size (tracked in CI artifacts, not a runtime benchmark)
- Cross-platform timing comparison (Linux vs macOS arm64)

**What is completely undefined:**

- Framework: `pytest-benchmark` is the obvious choice for Python-layer timing; `google/benchmark` for C++-layer microbenchmarks. Whether both are needed is unclear — if the Python layer is thin (which it is), Python-layer benchmarks probably suffice.
- CI vs. local: benchmark results are hardware-sensitive. Running in CI produces noisy numbers unless pinned to dedicated hardware. Options: (a) run locally on the owner's M-series Mac as the reference, with CI only checking that benchmarks don't error out; (b) use GitHub's larger runners for less-noisy CI results; (c) use a statistical approach (e.g., `pytest-benchmark`'s histogram mode) to tolerate noise.
- Acceptable regression thresholds: without a baseline, "10% regression" is undefined. Thresholds need to be set empirically after a first measurement pass.
- Whether benchmarks live in `tests/benchmarks/` (alongside pytest tests) or a separate `benchmarks/` directory.
- Whether benchmark results are committed to the repo (e.g., as JSON artifacts) or stored externally.

**Dependency chain:** Benchmarks are independent of items 1, 2, and 3. However, the most useful moment to establish a baseline is immediately before item 1 (OpenMesh → meshoptimizer), so the swap's performance impact can be quantified. Benchmarks should probably land in v0.4 as a prerequisite to item 1, not after.

---

## Dependencies

The items are not all independent:

```
Item 1: OpenMesh → meshoptimizer
    │
    └── unlocks Item 2: Meson migration
                (only worth evaluating if item 1 is done
                 AND CMake is still painful afterward)

Item 3: MCP server
    (independent — can proceed after Python API stabilizes in v0.3)

Item 4: Benchmarks
    (independent — ideally lands before item 1 to measure the decimation swap)
```

Recommended sequencing if all items are pursued:

1. Benchmarks (v0.4 or early v0.5) — establish baseline before architectural changes
2. OpenMesh → meshoptimizer (v0.5) — requires its own design pass; benchmarks confirm no regression
3. Meson migration (v0.5, maybe) — evaluate only after item 2; may be skipped entirely
4. MCP server (v0.5 or later) — independent, but depends on Python API stability

---

## ADR Implications

**D5 (FIRM — CMake + scikit-build-core):** The "what would invalidate" condition (OpenMesh dropped AND CMake pain measured) is the gate for item 2. D5 does not become FLEXIBLE until both conditions are met. If item 1 is done and CMake is fine, D5 stays FIRM and item 2 is dropped.

**D8 (FLEXIBLE — Keep OpenMesh for v0.1):** Item 1 is the execution of D8's deferred migration. D8 explicitly says it "deserves its own design pass" — that design pass is a prerequisite to any v0.5 implementation work. D8 does not become FIRM; it remains FLEXIBLE and its resolution is item 1's design doc.

**D1 (FIRM — C++20):** Unaffected. The half-edge-lite module for item 1 is C++20.

**D2 (FIRM — nanobind):** Unaffected by any of these items.

**D6 (FIRM — params/result structs):** Unaffected. The MCP server (item 3) would wrap the existing Python API, not redefine it.

**D7 (FLEXIBLE — harness):** Item 4 (benchmarks) extends the harness with a performance tier. The tier structure (Tier 1–5) defined in D7 covers correctness; benchmarks add a Tier 0 (performance baseline) or a separate performance harness. This needs to be specified in item 4's design pass.

---

## What's NOT in v0.5

Explicitly still not planned, regardless of what ships in v0.5:

- **GPU/CUDA backend** — no timeline. D1's rationale holds: MeshMonk is a batch CPU scientific tool. GPU registration adds substantial complexity with no clear demand from the user base.
- **Algorithmic improvements to registration math** — still "never, unless a measurable regression demands it" from the parent design doc. The half-edge-lite module in item 1 is implementation infrastructure, not an algorithmic change.
- **Rust rewrite** — D1 is FIRM. The conditions that would invalidate D1 (production-quality Rust halfedge library, full numerical re-derivation) have not been met.
- **Python API breaking changes** — v0.5 is an architectural layer below the Python API. The params/result struct API from D6 should survive intact.
- **MATLAB support** — D3 is FIRM. The university fork is the MATLAB path.
- **Windows-specific work beyond what v0.3 already ships** — v0.3 handles Windows wheels. v0.5 has no Windows-specific items.
- **Image-based visual regression tests** — still not planned. The Tier 3 human sign-off from v0.1 and the Tier 3.5 legacy-baseline gate remain the visual regression story.
