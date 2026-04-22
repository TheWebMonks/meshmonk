# API Polish Design

**Status:** Accepted  
**Date:** 2026-04-22  
**Decisions:** [ADR-003](../docs/decisions/ADR-003-api-polish.md)  
**Beads:** `meshmonk-modernization-hnf`, `-xvu`, `-h31`

## What We're Building

Three targeted improvements to the public Python API in `meshmonk/__init__.py`:

1. **hnf** — Consolidate duplicated Pattern A/B validation into `_prepare_arrays`
2. **xvu** — Standardise docstrings with a full kwarg→params mapping table
3. **h31** — Add complete type hints using `TypedDict` + `Unpack` (via `typing_extensions`)

All three are independent and can run in parallel. Each touches a different slice of `meshmonk/__init__.py` — the only coupling is merge-conflict risk in that file, resolvable by serialising merges.

## hnf: Validation Consolidation

Each of the three register functions currently duplicates an ~18-line guard that detects the call pattern and rejects mixed Pattern A/B usage. `_prepare_arrays` already owns pattern detection (`pattern_a = floating is not None`) and handles both paths — it's the right place to own validation too.

**What moves:** The Pattern A/B mixing guard and the `(floating is None) != (target is None)` symmetry check move into `_prepare_arrays`. Register functions drop those guard blocks and call `_prepare_arrays` directly.

**What stays:**
- Unknown-kwarg rejection stays in `_apply_rigid_kwargs` / `_apply_nonrigid_kwargs` / `_apply_pyramid_kwargs`. Those helpers raise `TypeError("rigid_register() got unexpected keyword arguments: ...")` with the function-named prefix — that prefix is user-facing and only the per-type helper knows the function name.
- `_mesh_to_arrays` continues to own the `normals_override` shape check. `_prepare_arrays` is not expanded to validate Pattern-A-only fields.

**Error-string contract:** The consolidated validator must preserve these exact strings (encoded as contracts in `tests/test_input_validation.py`):
- `"Cannot mix Pattern A (floating/target) with Pattern B (floating_features/target_features) arguments"` (6 tests match on `"Cannot mix Pattern A"`)
- `"Pattern B requires all array arguments. Missing: ..."` (2 tests match on `"Pattern B requires"`)
- `"both floating and target must be provided, or neither"` (symmetry check)

## xvu: Docstring Standardisation

`rigid_register` is the canonical reference function. Its docstring gets:
- Full Pattern A/B call pattern explanation (already present, keep)
- A Named-parameters subsection (explicit signature fields)
- A kwarg→params mapping table with per-function scope columns, replacing the placeholder `"See kwarg flattening rule in design doc"`

### Named parameters

These are explicit signature fields, not `**kwargs`. They appear in each register function's signature:

| name | type | scope |
|---|---|---|
| `floating`, `target` | `Mesh \| None` | Pattern A, all |
| `floating_features`, `target_features` | `NDArray \| None` | Pattern B, all |
| `floating_faces`, `target_faces` | `NDArray \| None` | Pattern B, all |
| `floating_flags`, `target_flags` | `NDArray \| None` | Pattern B, all |
| `normals`, `compute_normals_flag` | — | Pattern B only |
| `rigid_params` | `dict \| None` | nonrigid, pyramid |

### `**kwargs` parameters

Per-function scope uses Y/N columns instead of a verbose string, so the scope is unmistakable at a glance.

| kwarg | params field | rigid | nonrigid | pyramid |
|---|---|:-:|:-:|:-:|
| `correspondences_symmetric` | `params.correspondences.symmetric` | Y | Y | Y |
| `correspondences_num_neighbours` | `params.correspondences.num_neighbours` | Y | Y | Y |
| `correspondences_flag_threshold` ¹ | `params.correspondences.flag_threshold` | Y | Y | Y |
| `correspondences_equalize_push_pull` | `params.correspondences.equalize_push_pull` | Y | Y | Y |
| `inlier_kappa` | `params.inliers.kappa` | Y | Y | Y |
| `inlier_use_orientation` | `params.inliers.use_orientation` | Y | Y | Y |
| `num_iterations` | `params.num_iterations` | Y | Y | Y |
| `use_scaling` | `params.use_scaling` | Y | — | — |
| `transform_sigma` | `params.transform.sigma` | — | Y | Y |
| `transform_num_viscous_iterations_start` ² | `params.transform.num_viscous_iterations_start` | — | Y | Y |
| `transform_num_viscous_iterations_end` | `params.transform.num_viscous_iterations_end` | — | Y | Y |
| `transform_num_elastic_iterations_start` ² | `params.transform.num_elastic_iterations_start` | — | Y | Y |
| `transform_num_elastic_iterations_end` | `params.transform.num_elastic_iterations_end` | — | Y | Y |
| `downsample_float_start` | `params.downsample.float_start` | — | — | Y |
| `downsample_target_start` | `params.downsample.target_start` | — | — | Y |
| `downsample_float_end` | `params.downsample.float_end` | — | — | Y |
| `downsample_target_end` | `params.downsample.target_end` | — | — | Y |
| `num_pyramid_layers` | `params.num_pyramid_layers` | — | — | Y |

**Footnotes:**
1. `pyramid_register` defaults `correspondences_flag_threshold` to `0.999` (MATLAB convention) when the user doesn't pass it; `rigid_register` / `nonrigid_register` inherit the C++ struct default. The auto-population happens in `_apply_pyramid_kwargs` and relies on `kwargs.keys()` to detect explicit passes.
2. `pyramid_register` auto-populates `transform_num_viscous_iterations_start` and `transform_num_elastic_iterations_start` to `num_iterations` when not passed (MATLAB convention). `nonrigid_register` uses the C++ struct default.

`nonrigid_register` and `pyramid_register` add a single cross-reference line: *"Accepts the same two call patterns and `**kwargs` as `rigid_register`; see that docstring for the full parameter table."*

Passing a kwarg to a register function that doesn't accept it (e.g. `transform_sigma` to `rigid_register`) raises `TypeError` from the per-type `_apply_*_kwargs` helper; the table's Y/N columns are the source of truth for what's accepted.

## h31: Type Hints

Use `typing_extensions.Unpack` with `typing.TypedDict` to give `**kwargs` call-site type coverage without bumping the Python floor. This preserves the ADR-001 D2 commitment to CPython 3.10–3.13 and the 12-wheel abi3 matrix.

Define:

- `RigidKwargs(TypedDict, total=False)` — fields matching `_RIGID_KNOWN_KWARGS`
- `NonrigidKwargs(TypedDict, total=False)` — fields matching `_NONRIGID_KNOWN_KWARGS`
- `PyramidKwargs(TypedDict, total=False)` — fields matching `_PYRAMID_KNOWN_KWARGS`

Register function signatures use `**kwargs: Unpack[RigidKwargs]` etc. (imported from `typing_extensions`). This gives IDE autocomplete and pyright coverage at call sites on every supported Python version. `typing_extensions` is already a transitive dependency via `typer` and `pytest`; add it as a direct dependency.

**Runtime semantics to preserve:** `_apply_pyramid_kwargs` uses `kwargs.keys()` to detect explicitly-passed fields (vs defaults). `TypedDict(total=False)` preserves this — `.keys()` only enumerates present keys. Do not refactor to `kwargs.get(key, sentinel)` — that would collapse explicit-vs-default into a single branch and break pyramid's auto-population of `transform_num_viscous_iterations_start`, `transform_num_elastic_iterations_start`, and `correspondences_flag_threshold`.

**TypedDicts are not exported in `__all__`.** They live in a `meshmonk._types` submodule and are importable (`from meshmonk._types import RigidKwargs`) but not surfaced in the top-level public API. No downstream consumer has asked for them; keeping `__all__` lean leaves room to export later if a concrete need appears.

**Out of scope for h31:** Tightening result dataclass annotations from bare `np.ndarray` to `numpy.typing.NDArray[np.float32]` (and documenting shapes in docstrings). That's a separate, narrower task — filed as a follow-up bead, not rolled into h31.
