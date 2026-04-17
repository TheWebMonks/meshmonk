# Design: MeshMonk v0.2 — Team-ready

**Date:** 2026-04-17
**Status:** Draft
**Parent:** [Overview](./2026-04-17-meshmonk-modernization-design.md)
**Decisions:** [ADR-001](../docs/decisions/ADR-001-meshmonk-modernization.md)

---

## Scope

**v0.2 — Team-ready:** ship type stubs. README rewrite. Migration guide from MATLAB (side-by-side snippets) at `docs/migration-from-matlab.md`. Windows CI if feasible. cibuildwheel dry-run. Add diagnostic fields (`converged`, `fitness`, `inlier_rmse`, `num_inliers`) to result structs once convergence criteria are pinned down. PR `cli` → `master`; rename `master` → `main`.
