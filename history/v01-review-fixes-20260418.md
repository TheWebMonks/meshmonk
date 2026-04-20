# Fixes: v0.1 review findings implementation
Date: 2026-04-18
Review passes: 1 (architecture + implementation)

## Critical
(none)

## Important
- **library/src/meshmonk.cpp:~110** — Missing empty-faces guard in rigid_registration. nonrigid/pyramid have `if (floating_faces.rows() == 0 || target_faces.rows() == 0) return DegenerateInput;` but rigid does not. Fix: add the same check after line 111.

- **meshmonk/__init__.py:~498** — Silent flags loss in Pattern A. If user passes `floating_flags=custom` with Pattern A (`floating=mesh`), `floating_flags` is silently ignored. The mixed-pattern check only looks at `floating_features/target_features/floating_faces/target_faces`, not flags. Fix: add `floating_flags` and `target_flags` to the `pattern_b` detection, or warn when flags are provided with Pattern A.

- **library/src/NonrigidRegistration.cpp:~45** — NaN annealing rate when `num_viscous_iterations_start=0`. `exp(log(end/0)/(iters-1))` = Inf. PyramidNonrigidRegistration has a guard (`if start > 0`) but NonrigidRegistration does not. Fix: add the same guard.

- **tests/test_bug_fixes.py:291** — Hardcoded absolute path `/workspace/library/src/...`. Fix: use `Path(__file__).parent.parent / "library/src/PyramidNonrigidRegistration.cpp"`.

- **tests/test_api_enhancements.py:15** — Hardcoded Python version in LD_LIBRARY_PATH. Fix: compute dynamically via `sysconfig.get_paths()["platlib"]` or use existing `LD_LIBRARY_PATH` from environment.

## Minor
- **tests/test_api_enhancements.py:2** — Docstring says "v0.2" but this is v0.1 work.

- **tests/test_input_validation.py:8-15** — `_make_features` uses unseeded random. Fix: use deterministic vertices.

- **tests/test_registration.py:126-156** — `test_single_layer_matches_nonrigid` only checks both RMSEs < 1.0 independently, not that they're close to each other. Fix: add `assert abs(nr_rmse - pyr_rmse) < 0.5`.

- **tests/test_registration.py:181-186** — `test_translated_mesh_converges` only checks `aligned < original`. Fix: add magnitude check e.g. `assert aligned_dist < 0.5 * original_dist`.

## ADR Updates
- No ADR changes needed.

## Discarded
- **Duplicate Pattern B validation (~225 lines triplicated)**: Real observation but a v0.2 refactor, not a bug. Code works correctly.
- **MeshMonkError catches all std::runtime_error**: Known v0.1 constraint, ADR explicitly defers typed internal errors to v0.2.
- **std::cerr not routed through logger sink**: Known v0.1 constraint per ADR D6.
- **No column count validation on features before normals check**: Minor, C++ binding catches wrong-shape arrays independently.
- **demo --download stub exits code 0**: Pre-existing issue, not introduced by this change.
