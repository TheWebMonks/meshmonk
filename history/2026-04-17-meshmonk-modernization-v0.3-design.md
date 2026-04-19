# Design: MeshMonk v0.3 — PyPI-ready

**Date:** 2026-04-17
**Status:** Draft
**Parent:** [Overview](./2026-04-17-meshmonk-modernization-design.md)
**Decisions:** [ADR-001](../docs/decisions/ADR-001-meshmonk-modernization.md)

---

## Scope

**v0.3 — PyPI-ready:** cibuildwheel matrix (manylinux x86_64 + aarch64, macOS arm64, Windows). Stable ABI wheels for Python 3.12+ with per-version fallback wheels for 3.10/3.11. Release workflow with trusted-publisher. Test PyPI round-trip. MkDocs-material docs site. Publish to PyPI. `pip install meshmonk` works globally. Version bump to 0.3.0.

## Prerequisites

### PyPI trusted-publishing setup (manual, before first publish)

Before first publish: (a) claim the PyPI `meshmonk` name by uploading a placeholder `0.0.0.dev0` from `jsnyde0` manually; (b) configure trusted-publisher on `pypi.org/manage/account/publishing/` bound to `jsnyde0/meshmonk` + `.github/workflows/release.yml` + environment `pypi-release`; (c) dry-run on TestPyPI first. Workflow requires `permissions: id-token: write`.

### v0.2 deferred items disposition

Items explicitly deferred from v0.2 that v0.3 must address or defer further:

| Item | Disposition | Rationale |
|------|-------------|-----------|
| Diagnostic fields (`converged`, `fitness`, etc.) | **Defer to v0.4+** | Research/spec task — no convergence criterion exists in legacy code. Not blocking PyPI readiness. |
| `std::expected` migration | **Defer to v0.4+** | Raises compiler floor (gcc-13/clang-16), drops Ubuntu 22.04. `tl::expected` works correctly. Not blocking PyPI readiness. |
| Pluggable logger sink | **Defer to v0.4+** | `set_log_level("silent")` covers the primary use case. Pluggable output backend is a nice-to-have. |
| `RigidTransform.matrix` writable | **Defer to v0.4+** | Not blocking any v0.3 use case. |
| Validation boilerplate dedup (~69 lines) | **Defer to v0.4+** | Code quality, not functional. |
| `demo --download` actual URLs | **Include in v0.3** | v0.3 creates the first GitHub Release — natural place to wire up download. |
| Demo data redistribution rights | **Prerequisite check** | `demoFace.obj` has unknown license from "3Q Technologies Ltd." Demo data excluded from wheel, but GitHub Release hosting is still redistribution. Must resolve before hosting assets. If unresolved, `demo --download` ships as stub. |
| CI stub drift detection | **Defer to v0.4+** | Nice-to-have CI improvement. |
| `displacement_field` docstring | **Include in v0.3** | One-line fix, no risk. |

### ADR-001 D2 amendment required

ADR-001 D2 states "one wheel per OS/arch covers all supported Pythons (6 wheels total for v0.3: manylinux x86_64, manylinux aarch64, musllinux x86_64, macOS arm64, macOS x86_64, Windows x86_64)." The v0.3 design changes this to 12 wheels (4 abi3 + 8 per-version) and drops musllinux + macOS x86_64. The ADR must be amended before implementation starts:

- abi3 floor is Python 3.12 (not 3.10) due to nanobind's `PyType_FromMetaclass()` dependency
- Per-version wheels for Python 3.10/3.11 cover the gap
- musllinux dropped (niche audience, musl math-library differences)
- macOS x86_64 dropped (GitHub removed macOS 13 runners, Apple EOL)
- Total: 12 wheels (4 platforms × 3 Python versions, with abi3 for 3.12+)

---

## Stable ABI enablement

nanobind's stable ABI mode (`STABLE_ABI` flag) produces abi3 wheels. However, nanobind's STABLE_ABI requires **Python >= 3.12** because it depends on `PyType_FromMetaclass()` and limited-API vectorcall bindings added in CPython 3.12. On Python 3.10/3.11, the STABLE_ABI flag is silently ignored by nanobind's CMake, producing a normal (non-abi3) extension.

**Strategy:** Build `cp312-abi3` wheels (covering Python 3.12+) plus separate `cp310` and `cp311` per-version wheels for each platform. This gives full Python 3.10-3.13 coverage while using abi3 where possible. Python 3.13 is not listed in the build selector because the cp312-abi3 wheel is forward-compatible with 3.13+.

### CMake change

```cmake
# Before (current):
nanobind_add_module(_meshmonk_core bindings/bindings.cpp)

# After:
nanobind_add_module(_meshmonk_core STABLE_ABI bindings/bindings.cpp)
```

Note: nanobind automatically detects the Python version at build time. When building against Python 3.10/3.11, the STABLE_ABI flag is silently ignored and a normal extension is produced. When building against Python 3.12+, an abi3 extension is produced. This means the same CMakeLists.txt works for both abi3 and non-abi3 builds.

### pyproject.toml change

Do **NOT** add `py-api` to the static `[tool.scikit-build.wheel]` section — it would unconditionally tag all wheels (including 3.10/3.11 builds) as abi3. Instead, use a cibuildwheel override to set `py-api` only for 3.12 builds:

```toml
[[tool.cibuildwheel.overrides]]
select = "cp312-*"
environment = { SKBUILD_WHEEL_PY_API = "cp312" }
```

This ensures cp312 builds produce `cp312-abi3-*` tagged wheels, while cp310/cp311 builds produce standard per-version wheels.

### Verification

- Build abi3 wheel with Python 3.12, install on 3.12 and 3.13 — both must work
- Build per-version wheel with Python 3.10, install on 3.10 — must work
- Verify `nanobind_add_stub` still generates `.pyi` with `STABLE_ABI`
- Check wheel filenames: `*-cp312-abi3-*.whl` for 3.12+, `*-cp310-*.whl` for 3.10

### Known risks

If abi3 overhead exceeds ~2% (import time or registration latency), fall back to per-version wheels for all Pythons — no API impact, just remove the cibuildwheel override and STABLE_ABI flag.

---

## cibuildwheel matrix

### Target wheels (12 total: 4 abi3 + 8 per-version)

**abi3 wheels (Python 3.12+, 4 platforms):**

| OS | Arch | Selector | Notes |
|----|------|----------|-------|
| manylinux | x86_64 | `cp312-manylinux_x86_64` | `manylinux_2_28` (AlmaLinux 8) |
| manylinux | aarch64 | `cp312-manylinux_aarch64` | QEMU via `docker/setup-qemu-action`, `manylinux_2_28` |
| macOS | arm64 | `cp312-macosx_arm64` | macOS 14 runner |
| Windows | x86_64 | `cp312-win_amd64` | MSVC |

**Per-version wheels (Python 3.10, 3.11, 4 platforms each = 8 wheels):**

| OS | Arch | Selectors |
|----|------|-----------|
| manylinux | x86_64 | `cp310-manylinux_x86_64`, `cp311-manylinux_x86_64` |
| manylinux | aarch64 | `cp310-manylinux_aarch64`, `cp311-manylinux_aarch64` |
| macOS | arm64 | `cp310-macosx_arm64`, `cp311-macosx_arm64` |
| Windows | x86_64 | `cp310-win_amd64`, `cp311-win_amd64` |

**macOS x86_64 dropped:** GitHub removed macOS 13 runners (December 2025). Apple discontinued Intel Macs. macOS x86_64 users can build from source. Defer if demand materializes.

**musllinux dropped:** Niche audience for a 3D mesh registration library. musl has known math library differences. Defer to v0.4+ if demand materializes.

### CIBW configuration (pyproject.toml)

```toml
[tool.cibuildwheel]
build = "cp310-* cp311-* cp312-*"
test-requires = ["pytest", "scipy", "trimesh"]
test-command = "pytest {project}/tests -q --ignore={project}/tests/test_golden.py --ignore={project}/tests/test_silent_logger.py"

[tool.cibuildwheel.linux]
archs = ["x86_64", "aarch64"]
manylinux-x86_64-image = "manylinux_2_28"
manylinux-aarch64-image = "manylinux_2_28"

[tool.cibuildwheel.macos]
archs = ["arm64"]

[tool.cibuildwheel.windows]
archs = ["AMD64"]

[[tool.cibuildwheel.overrides]]
select = "cp312-*"
environment = { SKBUILD_WHEEL_PY_API = "cp312" }
```

`test_silent_logger.py` excluded because it spawns subprocesses that may not find the installed package inside cibuildwheel's ephemeral test environment. `test_golden.py` excluded because golden data is not shipped in sdist.

### CI workflow changes

- The existing `ci.yml` cibuildwheel dry-run image must be updated from `manylinux2014` to `manylinux_2_28` to match the release configuration
- `ci.yml` retains a single-platform cibuildwheel smoke test for PR validation
- `ci.yml` macos-13 matrix entries must be removed (runners deprecated December 2025)

---

## OpenMP handling

**Approach: verify, not fix.**

Investigation shows the vendored OpenMesh 11.0.0's `vci_openmp()` macro is defined but never invoked in the OpenMesh build tree. The `USE_OPENMP` cache variable is not checked by OpenMesh's CMake. OpenMP is likely not being linked.

### Verification steps

1. Build a wheel and check for OpenMP libraries:
   - Linux: `ldd _meshmonk_core*.so | grep omp`
   - macOS: `otool -L _meshmonk_core*.so | grep omp`
2. If no OpenMP linkage found: no action needed, document result
3. If OpenMP IS linked: investigate the actual mechanism and add the appropriate CMake fix

---

## Release workflow

New file: `.github/workflows/release.yml`

### Security: pin all actions to SHA

All third-party actions pinned to exact commit SHAs (not floating tags) to prevent supply-chain attacks. Comments note the version for readability:

```yaml
# Example (resolve actual SHAs at implementation time):
- uses: pypa/cibuildwheel@<sha>  # v2.x.y
- uses: pypa/gh-action-pypi-publish@<sha>  # release/v1
```

### Trigger

```yaml
on:
  push:
    tags: ["v*"]
  workflow_dispatch:
    inputs:
      publish_target:
        type: choice
        options: [testpypi, pypi, none]
        default: none
```

Tag push (`v0.3.0`, etc.) triggers full build + PyPI publish. Manual dispatch allows TestPyPI uploads or dry-run builds.

### Jobs

1. **build-wheels** — matrix of cibuildwheel targets (see matrix above). Each uploads wheel as artifact.
2. **build-sdist** — `python -m build --sdist`. Upload as artifact. Verify sdist contents match parent design spec (includes `library/`, `vendor/`, `CMakeLists.txt`, `pyproject.toml`; excludes `tests/`, `data/`, `history/`, `docs/`).
3. **test-sdist** — download sdist, `pip install` from source on Ubuntu with gcc, run tests. Verifies `pip install meshmonk --no-binary :all:` works. Note: Windows sdist build (`--no-binary :all:` with MSVC) is untested/best-effort.
4. **publish-testpypi** — (on `workflow_dispatch` with `testpypi`, or on pre-release tags `v*rc*`/`v*dev*`). Uses `pypa/gh-action-pypi-publish@<sha>` with TestPyPI. Requires `testpypi-release` environment.
5. **publish-pypi** — (on release tags `v*` without `rc`/`dev`). Uses `pypa/gh-action-pypi-publish@<sha>` with trusted-publisher. Requires `pypi-release` environment.
6. **create-release** — (on release tags). Creates GitHub Release with wheels + sdist as assets.

### Environment protection

GitHub environment `pypi-release`:
- Required reviewers: `jsnyde0`
- Deployment branches: tags matching `v*`

GitHub environment `testpypi-release`:
- No required reviewers (faster iteration during testing)
- Deployment branches: tags matching `v*`

### Permissions

```yaml
permissions:
  id-token: write    # trusted-publisher
  contents: write    # GitHub Release creation
```

---

## TestPyPI validation

### Process (before first real PyPI publish)

1. Bump version to `0.3.0.dev1` in `pyproject.toml`
2. Trigger release workflow manually with `publish_target: testpypi`
3. Wait for all wheels + sdist to build and upload
4. On a clean machine / fresh venv:
   ```bash
   pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ meshmonk
   python -c "import meshmonk; print(meshmonk.__version__)"
   ```
5. Run smoke tests (rigid registration on icosphere)
6. Verify abi3 wheel works on Python 3.12 and 3.13
7. Verify per-version wheel works on Python 3.10
8. If successful, bump to `0.3.0` final and tag for real release

### TestPyPI trusted-publisher

Separate trusted-publisher configuration on test.pypi.org:
- Same repo + workflow, but `testpypi-release` environment

---

## Documentation site (MkDocs-material)

### Scope note

The docs site is **not a blocker** for the PyPI publish. It can ship after the first release if needed. Minimum viable: landing page + installation + API reference.

### Structure

```
docs/
├── index.md                    # Landing page + installation
├── quickstart.md               # 5-minute rigid registration example
├── api/
│   ├── index.md                # API overview
│   ├── rigid.md                # rigid_register, compute_rigid_transform
│   ├── nonrigid.md             # nonrigid_register
│   ├── pyramid.md              # pyramid_register
│   ├── types.md                # RigidTransform, Result types, Params types
│   └── errors.md               # MeshMonkError, RegistrationError
├── migration-from-matlab.md    # Already written in v0.2
├── cli.md                      # CLI reference
└── decisions/                  # ADR summaries (link to full ADRs)
```

### Configuration

New file: `mkdocs.yml` at repo root.

```yaml
site_name: MeshMonk
theme:
  name: material
  palette:
    primary: blue
plugins:
  - search
  - mkdocstrings:
      handlers:
        python:
          options:
            show_source: false
```

### Dependencies

Add to `pyproject.toml`:

```toml
[project.optional-dependencies]
docs = ["mkdocs-material", "mkdocstrings[python]"]
```

### Hosting

GitHub Pages via `gh-pages` branch. CI deploys on tag push (part of release workflow) or manually via `mkdocs gh-deploy`.

---

## Version management

### Single source of truth

`pyproject.toml` `version = "0.3.0"` is the single source. `meshmonk/__init__.py` reads it via:

```python
try:
    from importlib.metadata import version
    __version__ = version("meshmonk")
except Exception:
    __version__ = "0.0.0.dev0"  # fallback for uninstalled dev runs
```

This replaces the current hardcoded `__version__` string and eliminates the dual-maintenance problem. The `try/except` handles the case where tests run directly from the repo without `pip install -e .`.

### Bump timing

- First v0.3 commit: bump to `0.3.0.dev0`
- TestPyPI validation: `0.3.0.dev1`, `0.3.0rc1`, etc.
- Final release: `0.3.0`
- Tag: `v0.3.0`

---

## `demo --download` implementation

Since v0.3 creates the first GitHub Release, this is when `demo --download` can be wired up.

### Implementation

```python
# meshmonk/cli.py
DEMO_ASSETS = {
    "Template.obj": {
        "url": "https://github.com/jsnyde0/meshmonk/releases/latest/download/Template.obj",
        "sha256": "<computed-at-release-time>",
    },
    "demoFace.obj": {
        "url": "https://github.com/jsnyde0/meshmonk/releases/latest/download/demoFace.obj",
        "sha256": "<computed-at-release-time>",
    },
}
CACHE_DIR = Path("~/.cache/meshmonk").expanduser()
```

Uses `/releases/latest/download/` URLs so new releases automatically serve current assets without code changes. SHA-256 pins specific content for integrity verification. On SHA-256 mismatch, raise an error with message directing users to update their meshmonk version (safest behavior — prevents running with potentially incompatible demo data).

### Prerequisite

Demo data redistribution rights must be resolved first (see Prerequisites section). If `demoFace.obj` cannot be redistributed, use only `Template.obj` or provide a synthetic demo mesh.

---

## Implementation ordering

### Phase 1 — Build infrastructure (no functional changes)

1. Update ci.yml: remove macos-13 entries (runners deprecated December 2025), update cibuildwheel dry-run image from `manylinux2014` to `manylinux_2_28`
2. Amend ADR-001 D2 with updated wheel strategy (12 wheels, abi3 floor cp312, dropped platforms)
3. Stable ABI enablement (CMakeLists.txt `STABLE_ABI` flag)
4. pyproject.toml updates (cibuildwheel config with cp312 override, docs extra)
5. OpenMP verification (check wheel contents, document result)
6. Version management migration (`importlib.metadata` with fallback)

### Phase 2 — Release workflow

7. Create `.github/workflows/release.yml` with SHA-pinned actions
8. TestPyPI trusted-publisher setup (manual)
9. TestPyPI round-trip validation
10. PyPI trusted-publisher setup (manual)

### Phase 3 — Documentation

11. MkDocs site setup (`mkdocs.yml`, initial pages)
12. API reference generation
13. GitHub Pages deployment in release workflow

### Phase 4 — Release

14. `demo --download` wiring (if redistribution rights resolved)
15. `displacement_field` docstring fix
16. Version bump to `0.3.0` final
17. Tag and publish

### Manual steps (owner required)

- PyPI name claim (placeholder upload)
- Trusted-publisher configuration on pypi.org and test.pypi.org
- GitHub environment creation (`pypi-release`, `testpypi-release`)
- Demo data redistribution rights resolution
- Final visual check of docs site
