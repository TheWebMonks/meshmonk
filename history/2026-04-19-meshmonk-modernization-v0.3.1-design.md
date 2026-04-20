# Design: MeshMonk v0.3.1 — Ship It

**Date:** 2026-04-19
**Status:** Draft
**Parent:** [Overview](./2026-04-17-meshmonk-modernization-design.md)
**Decisions:** [ADR-001](../docs/decisions/ADR-001-meshmonk-modernization.md)

---

## Scope

Pre-merge housekeeping before (or alongside) the PR of `meshmonk-modernization` → `master`. Six independent items: one branch operation, one repo hygiene change, one comment fix, and three CI additions. All items are small and well-understood. No new design is required; this doc is structured tracking only.

## Prerequisites

None of the items below require owner-level prerequisites (no PyPI name claim, no OIDC trusted-publisher config, no GitHub Environments setup, no demo data rights, no TestPyPI round-trip). The branch rename (item 2) requires the repository owner to update the GitHub default-branch setting in repo Settings after the git-level rename.

## Items

### 1. PR meshmonk-modernization → master

The v0.2 design doc (line 12) planned: PR `meshmonk-modernization` → `cli` → `master`. The `cli` intermediate step is skippable. The `cli/` directory was deleted as dead code (v0.2 Priority 1, item 3), and no separate `cli` branch carries work that is not already on `meshmonk-modernization`. The 133-commit branch can go directly to `master`.

**Action:** Open a PR from `meshmonk-modernization` to `master`. Confirm `ci.yml` triggers on the PR (it already lists `master` in `pull_request.branches`).

### 2. Rename master → main

Planned in v0.2 design doc (line 12: "rename `master` → `main`"), never done. Standard modern convention; GitHub now defaults new repos to `main`.

**Action:**
1. After the PR in item 1 merges, rename the branch:
   ```
   git branch -m master main
   git push origin -u main
   git push origin --delete master
   ```
2. In GitHub repo Settings → Branches, update the default branch to `main`.
3. Update `ci.yml` `pull_request.branches` — currently lists `["master", "main"]`, so `main` is already covered. Remove the now-dead `"master"` entry.
4. Check `release.yml` for any hardcoded references to `master` (none observed, but verify).

**Note:** Items 3–6 should target `main` after this rename. If done before the rename, they land on `master` (which becomes `main`), so ordering does not affect the final state.

### 3. Fix stale TODO comments

Three comments in `library/src/meshmonk.cpp` (lines 165, 257, 366) read:

```cpp
// TODO (v0.3): add convergence criterion (NonConvergence error)
```

Convergence was explicitly deferred past v0.3. The v0.2 fixes file (line 26) updated `TODO (v0.2)` → `TODO (v0.3)` for these same comments in an earlier round, but the v0.3 design explicitly deferred convergence again. The correct label is now `v0.4`.

**Action:** Change all three instances to:
```cpp
// TODO (v0.4): add convergence criterion (NonConvergence error)
```

This is a one-liner per site; no tests required (comment-only change).

### 4. GH Pages deployment

The `docs/` directory and `mkdocs.yml` exist and the mkdocs site builds locally. No CI job deploys it. The v0.3 fixes file (discarded section) noted: "design doc explicitly says docs site 'can ship after the first release.'" That time has come — the first release is imminent.

**Action:** Add a `deploy-docs` job to `ci.yml` (or a standalone `docs.yml` workflow) that runs `mkdocs gh-deploy --force` on push to `main`. Use the standard `peaceiris/actions-gh-pages` or the built-in `mkdocs gh-deploy` approach with `actions/checkout` + `mkdocs build`. The job should only run on `push` to `main`, not on PRs.

**SHA-pinning:** Since item 6 will SHA-pin `ci.yml` actions anyway, coordinate — pin the new docs action at the same time.

### 5. mypy/pyright CI

The v0.2 design specified a type-checking CI job. The v0.2 review fixes file (line 55) noted: "Design doc specifies mypy/pyright CI, not implemented" and deferred it. The project already ships type stubs (`meshmonk/_meshmonk_core.pyi`) and a `py.typed` marker; type checking is now meaningful.

**Action:** Add a `typecheck` job to `ci.yml`. Suggested form:

```yaml
typecheck:
  name: typecheck (ubuntu-latest, python 3.12)
  runs-on: ubuntu-latest
  steps:
    - uses: actions/checkout@<SHA>  # see item 6
    - uses: actions/setup-python@<SHA>
      with:
        python-version: "3.12"
    - name: Install uv
      run: pip install uv
    - name: Install package with dev extra
      run: uv pip install --system ".[dev]"
    - name: Run pyright
      run: pyright meshmonk/
```

Confirm `pyright` is listed in `[project.optional-dependencies] dev` in `pyproject.toml`; add it if missing.

### 6. ci.yml SHA-pinning

`release.yml` uses SHA-pinned actions throughout (e.g., `actions/checkout@11bd71901bbe5b1630ceea73d27597364c9af683  # v4.2.2`). `ci.yml` uses floating version tags (`@v4`, `@v5`, `@v3.28.0`, `@v1`, `@v2`). The v0.3 fixes file (discarded section) noted this gap explicitly. Supply-chain consistency requires both files to pin.

**Current unpinned actions in ci.yml:**

| Action | Current ref | Pinned equivalent (resolve before committing) |
|--------|-------------|-----------------------------------------------|
| `actions/checkout` | `@v4` | SHA for v4.2.2 (same as release.yml: `11bd71901bbe5b1630ceea73d27597364c9af683`) |
| `actions/setup-python` | `@v5` | SHA for v5.3.0 (same as release.yml: `0b93645e9fea7318ecaed2b359559ac225c90a2b`) |
| `lukka/get-cmake` | `@v3.28.0` | Resolve SHA for v3.28.0 via `git ls-remote https://github.com/lukka/get-cmake refs/tags/v3.28.0` |
| `ilammy/msvc-dev-cmd` | `@v1` | Resolve SHA for latest v1.x tag |
| `pypa/cibuildwheel` | `@v2` | SHA for v2.22.0 (same as release.yml: `ee63bf16da6cddfb925f542f2c7b59ad50e93969`) |
| `actions/upload-artifact` | `@v4` | SHA for v4.5.0 (same as release.yml: `6f51ac03b9356f520e9adb1b1b7802705f340c2b`) |

**Action:** For each action, resolve the exact commit SHA for the version tag in use and replace `@vX[.Y.Z]` with `@<SHA>  # vX.Y.Z`. Reuse the SHAs already present in `release.yml` where the version matches.

## Ordering

All six items are independent of each other with one exception:

- **Item 1 must complete before item 2** (can't rename to `main` before the merge-to-master PR lands, or the PR target disappears).
- Items 3, 4, 5, 6 are fully independent of each other and of item 1. They can be done before or after the PR merges; they will land on `main` either way.
- Items 4 and 6 share a coordination point: if both touch `ci.yml` in the same commit, do them together to avoid a redundant pin-then-add churn.

Suggested sequence:

1. Items 3, 5, 6 (and 4 if coordinating with 6) — land on `meshmonk-modernization` now
2. Item 1 — open and merge the PR
3. Item 2 — rename branch after merge

## Risks

- **Item 2 (branch rename):** Any clone with a cached `origin/master` remote reference will break until `git fetch --prune`. CI secrets, branch-protection rules, and any external integrations (webhooks, badges) that hardcode `master` will need updating. Low risk for a single-maintainer repo but requires a checklist pass.
- **Item 4 (GH Pages):** If the `gh-pages` branch does not yet exist, `mkdocs gh-deploy` creates it on first run. No destructive risk, but the first deploy may require the repo to have GitHub Pages enabled in Settings (Source: Deploy from branch, branch: `gh-pages`).
- **Item 5 (pyright CI):** The type stubs cover the C extension surface; the pure-Python layer in `meshmonk/__init__.py` may have latent type errors not previously caught. The job might fail on first run. Treat initial failures as findings to fix, not as reasons to skip the job.
- **Item 6 (SHA-pinning):** `lukka/get-cmake@v3.28.0` and `ilammy/msvc-dev-cmd@v1` are not already pinned in `release.yml`, so their SHAs must be resolved fresh. Use `git ls-remote` against the upstream repo to get the tag's commit SHA.
