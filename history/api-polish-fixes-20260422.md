# Fixes: api-polish
Date: 2026-04-22
Review passes: 1 (architecture + implementation in parallel)

Commits reviewed: fed7e9c (hnf), eb5f7ed (xvu), 649c619 (h31), 7abdf9f (epic-closure docs).

## Critical
*None.*

## Important

- **meshmonk/__init__.py:623-627 (rigid_register docstring) + history/2026-04-22-api-polish-design.md:50** — `normals` and `compute_normals_flag` are labeled "Pattern B only" but are actually **Pattern A only**. Both are only consumed inside the `if pattern_a:` branch of `_prepare_arrays` (lines 507-514 in __init__.py); the Pattern B branch discards them silently. A user calling with Pattern B arrays and `normals=my_normals` will get a silent no-op. Change the docstring subsection header from "Pattern B only:" to "Pattern A only:" and fix the same mislabel in the design doc.

- **meshmonk/__init__.py:629-631 (rigid_register docstring)** — `rigid_params` is documented under `rigid_register`'s Named Parameters section, but `rigid_register` does not accept it (it raises `TypeError: rigid_register() got unexpected keyword arguments: ['rigid_params']` via `_apply_rigid_kwargs`). Remove the entire "nonrigid and pyramid only: rigid_params" block from `rigid_register`'s docstring. `rigid_params` is already documented in `nonrigid_register` (lines 730-734) and `pyramid_register` (lines 816-821) where it belongs; the cross-reference line in those docstrings covers the shared `**kwargs` table.

- **pyproject.toml:16** — `typing_extensions>=4.0` is below PEP 692's `Unpack[TypedDict]` checker support. `Unpack` exists in typing_extensions from 4.1.0, basic `Unpack[TypedDict]` is 4.3.0, and pyright/mypy standardized recognition is 4.6.0. Bump to `typing_extensions>=4.6` to guarantee the call-site type coverage that motivated h31. This is not an ADR-003 D3 conflict — D3 prescribed the approach, not a specific floor.

- **meshmonk/__init__.py:421-422, 451-454, 866** — `_apply_pyramid_kwargs` accepts `explicit_kwargs: set[str]` as a third parameter, but the only production call site (line 866) always passes `set(kwargs.keys())`, and the test at `tests/test_kwargs_refactor.py:240` demonstrates the parameter can drift from truth (passes `explicit_kwargs=set()` while `kwargs={"bogus_param":1}`). Remove the `explicit_kwargs` parameter and replace `"... not in explicit_kwargs"` with `"... not in kwargs"` on lines 451 and 453 — matching the pattern already used one layer up on line 864. Update the 6 affected call sites in `tests/test_kwargs_refactor.py` (lines 80-83, 120-123, 168-175, 197, 210, 240) and the 2 in `tests/test_type_hints_h31.py`.

- **meshmonk/__init__.py:451-454, 861-865** — Pyramid's three MATLAB-convention "default if not explicitly passed" auto-populations are split across two helpers. `correspondences_flag_threshold` (→ 0.999) lives in `pyramid_register`'s body (line 861-865); `transform_num_viscous_iterations_start` and `transform_num_elastic_iterations_start` (→ `params.num_iterations`) live in `_apply_pyramid_kwargs` (lines 451-454). These are the same conceptual mechanism and share footnotes in the docstring table. Move the `correspondences_flag_threshold = 0.999` auto-populate into `_apply_pyramid_kwargs` alongside the other two, so all pyramid MATLAB defaults live in one function and `pyramid_register`'s body matches `rigid_register` / `nonrigid_register` in shape.

- **meshmonk/__init__.py:721, 804** — `rigid_params: dict | None = None` in `nonrigid_register` and `pyramid_register` signatures leaves an untyped hole in the otherwise TypedDict-annotated API. `RigidKwargs` already exists from h31. Change the parameter annotation to `rigid_params: RigidKwargs | None = None` in both function signatures and update the docstring types. The runtime `isinstance(rigid_params, dict)` guards (lines 756, 842) still work because TypedDicts are dicts at runtime. Also update the Named Parameters docstring entries in both functions to reference `RigidKwargs`.

- **tests/test_type_hints_h31.py:276-298** — `test_rigid_register_runtime_behavior_unchanged` wraps the full `rigid_register` call in `try: ... except Exception: pass`, so any crash, wrong kwarg handling, or silent no-op passes the test. The docstring claims it verifies "Kwargs still apply correctly at runtime" but it verifies nothing. Remove the test. If runtime kwarg-application coverage is wanted, `tests/test_kwargs_refactor.py` already exercises `_apply_rigid_kwargs` directly against the params struct — that's the right layer for this check.

## Minor

- **meshmonk/__init__.py:352, 370, 392, 407, 421-422** — The `Mapping[str, Any]` annotation on the five `_apply_*_kwargs` helper signatures is load-bearing (TypedDicts are subtypes of `Mapping[str, object]`, not `dict`). A future contributor "cleaning up" back to `dict` would silently break pyright at the TypedDict call sites. Add a one-line comment near `_apply_shared_kwargs` (or adjacent to the `Mapping` import) noting this: `# Mapping[str, Any] accepts both plain dict and TypedDict; do not narrow to dict.`

- **tests/test_type_hints_h31.py:322-383** — Field-type tests (`test_rigid_kwargs_field_types` and siblings) access `TypedDict.__annotations__` directly and assert identity against the bare type object (`annotations["correspondences_symmetric"] is bool`). This only works because `_types.py` lacks `from __future__ import annotations`; adding that import (a natural refactor since `__init__.py` already has it) would silently break all these tests because annotations would be strings. Replace direct `__annotations__` access with `typing.get_type_hints(TypedDictClass)` — which resolves string annotations to actual type objects regardless of the future-annotations setting. This matches the pattern already used in `test_rigid_register_kwargs_annotation` (line 424), which uses `get_type_hints(..., include_extras=True)`.

## ADR Updates

- **ADR-003 D3**: no update needed — the decision stands; the fix is tightening the version floor to reflect the already-stated rationale (call-site type coverage). `typing_extensions>=4.6` better captures D3's intent without changing the decision.

- **Design doc `history/2026-04-22-api-polish-design.md:50`**: fix the `normals` / `compute_normals_flag` mislabel as part of fix #1 above.

## Discarded

- **Architecture finding 7 (`meshmonk/_types.py` needs `__all__` or re-export from `meshmonk/__init__.py`):** The design doc at line 100 explicitly documents `from meshmonk._types import RigidKwargs` as the supported import path — "importable but not surfaced in the top-level public API." The private-submodule path is intentional, not an oversight. Re-exporting from `meshmonk` would contradict the documented design.

- **cp313 wheel matrix, pyright/mypy config absent, result dataclass dtype tightening (follow-up meshmonk-modernization-uij.1):** all verified as non-issues for this change. uij.1 is already filed as an intentional follow-up per the design doc's h31 scope decision.

- **`uv.lock` diff scope:** verified as minimal (typing-extensions addition only; no unrelated package upgrades). No fix needed.
