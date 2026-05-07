#!/usr/bin/env bash
# Run the e2e demo inside a Linux container for visual smoke testing.
#
# Two modes:
#   wheel  (default, ~30s) — install latest dev wheel from TestPyPI
#   src    (~3 min)        — build current working tree from source
#                            (use this to verify uncommitted/unreleased changes)
#
# Usage:
#   bash scripts/linux_e2e_smoke.sh         # wheel mode
#   bash scripts/linux_e2e_smoke.sh src     # source build
#
# Outputs PNGs + OBJs to ./demo/linux_out/ on the host.
set -euo pipefail

MODE="${1:-wheel}"
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT_DIR="$REPO_ROOT/demo/linux_out"
mkdir -p "$OUT_DIR"

echo ">>> Mode: $MODE"
echo ">>> Repo: $REPO_ROOT"
echo ">>> Output dir (host): $OUT_DIR"

case "$MODE" in
  wheel)
    # Fast path: pip-install the latest TestPyPI wheel, then run demo.
    # Demo data (Template.obj, demoFace.obj) lives inside demo/.
    docker run --rm \
      -v "$REPO_ROOT/demo":/demo:ro \
      -v "$OUT_DIR":/out \
      -w /work \
      python:3.12-bookworm \
      bash -euxc '
        cp -a /demo/. /work/

        pip install --quiet --upgrade pip
        pip install --quiet \
          --index-url https://test.pypi.org/simple/ \
          --extra-index-url https://pypi.org/simple/ \
          meshmonk trimesh matplotlib scipy

        python e2e_demo.py

        cp -v e2e_compare.png e2e_rigid.png e2e_nonrigid.png e2e_pyramid.png /out/
        cp -v e2e_rigid.obj  e2e_nonrigid.obj  e2e_pyramid.obj  /out/
      '
    ;;
  src)
    # Slow path: build current working tree from source. Use this when
    # validating uncommitted/unreleased changes.
    docker run --rm \
      -v "$REPO_ROOT":/src:ro \
      -v "$OUT_DIR":/out \
      -w /work \
      python:3.12-bookworm \
      bash -euxc '
        apt-get update -qq
        apt-get install -y -qq cmake build-essential >/dev/null

        cp -a /src/. /work/

        pip install --quiet --upgrade pip
        pip install --quiet ".[io]" matplotlib scipy

        cd demo
        python e2e_demo.py

        cp -v e2e_compare.png e2e_rigid.png e2e_nonrigid.png e2e_pyramid.png /out/
        cp -v e2e_rigid.obj  e2e_nonrigid.obj  e2e_pyramid.obj  /out/
      '
    ;;
  *)
    echo "Unknown mode: $MODE (expected 'wheel' or 'src')"
    exit 1
    ;;
esac

echo
echo ">>> Done. View PNGs:"
echo "    open $OUT_DIR/e2e_compare.png"
ls -la "$OUT_DIR"
