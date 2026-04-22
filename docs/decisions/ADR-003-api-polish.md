# ADR-003: API Polish — Validation Consolidation, Docstrings, Type Hints

**Status:** Accepted  
**Date:** 2026-04-22  
**Design:** [API polish design doc](../../history/2026-04-22-api-polish-design.md)  
**Parent:** [ADR-001 meshmonk-modernization](./ADR-001-meshmonk-modernization.md)  
**Related:** Beads `meshmonk-modernization-hnf`, `-xvu`, `-h31`

## Context

Three polish tasks were filed after a review of the public Python API in `meshmonk/__init__.py`. Each is narrow but the decisions ripple across the others: validation shape (hnf) determines what callers look like, docstrings (xvu) document that shape, and type hints (h31) formalise it.

## Decisions

### D1: Merge Pattern A/B validation into `_prepare_arrays`

**Firmness: FIRM**

The 18-line guard that detects and rejects mixed Pattern A/B usage is currently copied verbatim into all three register functions. `_prepare_arrays` already owns pattern detection and handles both paths — validation belongs there too.

**Rationale:** Single ownership. `_prepare_arrays` is the only function that sees the full argument set and therefore the only place that can make the determination cleanly. Callers should not know about patterns.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **Merge into `_prepare_arrays` (chosen)** | Single owner, callers simplified | `_prepare_arrays` slightly heavier |
| Extract as `_validate_call_pattern()` helper | Explicit, testable in isolation | Two callsites per register function instead of one |
| Leave in callers | No change | Drift risk as new register functions are added |

**What would invalidate this:** A future register function that needs different validation semantics — at that point a separate helper may be warranted.

---

### D2: `rigid_register` is the docstring source of truth

**Firmness: FLEXIBLE**

The full Pattern A/B explanation and kwarg→params mapping table live in `rigid_register`. `nonrigid_register` and `pyramid_register` cross-reference it with one line.

**Rationale:** `rigid_register` is the simplest and most commonly read function. Full duplication across three docstrings means three places to update when a kwarg is added. Templating machinery solves the same problem with more complexity than the problem warrants.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **Cross-reference pattern (chosen)** | Zero machinery, easy to maintain | Caller must navigate to `rigid_register` for full docs |
| Shared string constant | Docs stay in sync automatically | String interpolation obscures docstrings in IDEs |
| Programmatic assembly | Same as above | Same cons, more code |

**What would invalidate this:** If `nonrigid_register` or `pyramid_register` develop meaningfully different call semantics that make cross-referencing misleading.

---

### D3: Use `typing_extensions.Unpack` + `TypedDict` for kwargs, keep 3.10 floor

**Firmness: FLEXIBLE**

Define `RigidKwargs`, `NonrigidKwargs`, `PyramidKwargs` as `TypedDict`s and use `Unpack` from `typing_extensions` in register function signatures. Keep `requires-python = ">=3.10"`. TypedDicts live in `meshmonk._types` and are importable but not added to `__all__`.

**Rationale:** `typing_extensions.Unpack` backports PEP 692 to Python 3.10/3.11 with full pyright and mypy support, so the call-site UX win is available without bumping the Python floor. `typing_extensions` is already a transitive dependency (via `typer`, `pytest`); promoting it to a direct dep is a one-line change. Keeping the floor at 3.10 preserves ADR-001 D2 (FIRM) — the 12-wheel CPython 3.10–3.13 matrix — and avoids dropping Ubuntu 22.04 LTS (3.10) and Debian 12 (3.11) users on a 0.3.0.dev library.

**Alternatives considered:**

| Approach | Pros | Cons |
|---|---|---|
| **`typing_extensions.Unpack`, 3.10 floor (chosen)** | Full IDE/pyright coverage at call sites; preserves ADR-001 D2; no user-visible breaking change | One extra direct dep (already transitively present) |
| Bump floor to 3.12, use `typing.Unpack` | Stdlib only | CONFLICTS WITH ADR-001 D2 (FIRM); drops cp310/cp311 wheels; collapses 12-wheel matrix to 4; visible compat break on pre-1.0 library |
| Stay on 3.10, `**kwargs: Any` | Broader compat, zero new deps | No call-site type benefit |
| Stay on 3.10, use `@overload` stubs | Call-site coverage without new dep | Significant boilerplate, fragile to maintain |

**What would invalidate this:** `typing_extensions` dropping `Unpack` support, or ADR-001 D2 being amended to drop 3.10/3.11 for other reasons — at which point we can switch to `typing.Unpack`.
