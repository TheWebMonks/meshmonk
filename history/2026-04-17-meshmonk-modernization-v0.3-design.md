# Design: MeshMonk v0.3 — PyPI-ready

**Date:** 2026-04-17
**Status:** Draft
**Parent:** [Overview](./2026-04-17-meshmonk-modernization-design.md)
**Decisions:** [ADR-001](../docs/decisions/ADR-001-meshmonk-modernization.md)

---

## Scope

**v0.3 — PyPI-ready:** cibuildwheel matrix (manylinux, macOS x86_64 + arm64, Windows). Test PyPI round-trip. MkDocs-material docs site. Publish to PyPI. `pip install meshmonk` works globally.

## PyPI trusted-publishing setup

- **PyPI trusted-publishing setup (v0.3 prerequisite):** before first publish: (a) claim the PyPI `meshmonk` name by uploading a placeholder `0.0.0.dev0` from `jsnyde0` manually; (b) configure trusted-publisher on `pypi.org/manage/account/publishing/` bound to `jsnyde0/meshmonk` + `.github/workflows/release.yml` + environment `pypi-release`; (c) dry-run on TestPyPI first. Workflow requires `permissions: id-token: write`.
