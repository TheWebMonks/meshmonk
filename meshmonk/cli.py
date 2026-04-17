"""MeshMonk command-line interface (typer-based).

Entry point: meshmonk.cli:app  (registered in pyproject.toml)

Subcommands
-----------
rigid      -- run rigid (SE(3)) registration on two OBJ meshes
nonrigid   -- run nonrigid (viscoelastic) registration
pyramid    -- run pyramid (multi-resolution nonrigid) registration
demo       -- download demo meshes and/or run a demo registration

All mesh I/O is handled by trimesh (meshmonk[io] extra).  If trimesh is not
installed the CLI emits a clear error message rather than a raw traceback.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

app = typer.Typer(
    name="meshmonk",
    help="MeshMonk 3D mesh registration CLI.",
    no_args_is_help=True,
)

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_TRIMESH_IMPORT_ERROR = (
    "trimesh is required for mesh I/O.  Install it with:\n"
    "  pip install 'meshmonk[io]'"
)

# Default download URLs — set to None as a placeholder until hosting is confirmed.
_TEMPLATE_URL: Optional[str] = None
_DEMO_FACE_URL: Optional[str] = None

_CACHE_DIR = Path.home() / ".cache" / "meshmonk"
_REPO_DATA_DIR = Path(__file__).parent.parent / "data"


def _require_trimesh():
    """Import and return trimesh, or exit with a helpful message."""
    try:
        import trimesh  # noqa: PLC0415
        return trimesh
    except ImportError:
        typer.echo(_TRIMESH_IMPORT_ERROR, err=True)
        raise typer.Exit(code=1)


def _load_mesh(trimesh_mod, path: str):
    """Load an OBJ (or any trimesh-supported format) from *path*."""
    return trimesh_mod.load(str(path))


def _save_mesh(trimesh_mod, vertices, faces, out_path: str) -> None:
    """Save a mesh defined by *vertices* + *faces* to *out_path*."""
    import numpy as np  # noqa: PLC0415
    mesh = trimesh_mod.Trimesh(
        vertices=np.asarray(vertices, dtype=np.float64),
        faces=np.asarray(faces, dtype=np.int32),
        process=False,
    )
    mesh.export(str(out_path))


def _find_demo_meshes():
    """Return (template_path, demo_face_path) or None if not found.

    Search order:
    1. ~/.cache/meshmonk/ (downloaded via ``meshmonk demo --download``)
    2. <repo-root>/data/ (developer convenience)
    """
    for search_dir in (_CACHE_DIR, _REPO_DATA_DIR):
        t = search_dir / "Template.obj"
        d = search_dir / "demoFace.obj"
        if t.exists() and d.exists():
            return t, d
    return None


# ---------------------------------------------------------------------------
# rigid subcommand
# ---------------------------------------------------------------------------

@app.command()
def rigid(
    floating: str = typer.Argument(..., help="Path to floating (source) mesh OBJ."),
    target: str = typer.Argument(..., help="Path to target mesh OBJ."),
    out: str = typer.Option("result.obj", "--out", help="Output OBJ path."),
    iterations: int = typer.Option(80, "--iterations", help="Number of ICP iterations."),
    use_scaling: bool = typer.Option(False, "--use-scaling/--no-use-scaling", help="Allow uniform scaling."),
    kappa: float = typer.Option(12.0, "--kappa", help="Inlier detection kappa."),
    num_neighbours: int = typer.Option(3, "--num-neighbours", help="Correspondence neighbours."),
) -> None:
    """Run rigid (SE(3)) mesh registration on FLOATING → TARGET."""
    trimesh = _require_trimesh()
    import meshmonk  # noqa: PLC0415

    floating_mesh = _load_mesh(trimesh, floating)
    target_mesh = _load_mesh(trimesh, target)

    result = meshmonk.rigid_register(
        floating=floating_mesh,
        target=target_mesh,
        num_iterations=iterations,
        use_scaling=use_scaling,
        inlier_kappa=kappa,
        correspondences_num_neighbours=num_neighbours,
    )

    _save_mesh(trimesh, result.aligned_vertices, floating_mesh.faces, out)
    typer.echo(f"Saved result to {out}")


# ---------------------------------------------------------------------------
# nonrigid subcommand
# ---------------------------------------------------------------------------

@app.command()
def nonrigid(
    floating: str = typer.Argument(..., help="Path to floating (source) mesh OBJ."),
    target: str = typer.Argument(..., help="Path to target mesh OBJ."),
    out: str = typer.Option("result.obj", "--out", help="Output OBJ path."),
    iterations: int = typer.Option(200, "--iterations", help="Number of iterations."),
    kappa: float = typer.Option(12.0, "--kappa", help="Inlier detection kappa."),
    num_neighbours: int = typer.Option(3, "--num-neighbours", help="Correspondence neighbours."),
) -> None:
    """Run nonrigid (viscoelastic) mesh registration on FLOATING → TARGET."""
    trimesh = _require_trimesh()
    import meshmonk  # noqa: PLC0415

    floating_mesh = _load_mesh(trimesh, floating)
    target_mesh = _load_mesh(trimesh, target)

    result = meshmonk.nonrigid_register(
        floating=floating_mesh,
        target=target_mesh,
        num_iterations=iterations,
        inlier_kappa=kappa,
        correspondences_num_neighbours=num_neighbours,
    )

    _save_mesh(trimesh, result.aligned_vertices, floating_mesh.faces, out)
    typer.echo(f"Saved result to {out}")


# ---------------------------------------------------------------------------
# pyramid subcommand
# ---------------------------------------------------------------------------

@app.command()
def pyramid(
    floating: str = typer.Argument(..., help="Path to floating (source) mesh OBJ."),
    target: str = typer.Argument(..., help="Path to target mesh OBJ."),
    out: str = typer.Option("result.obj", "--out", help="Output OBJ path."),
    iterations: int = typer.Option(90, "--iterations", help="Total iterations across all pyramid layers."),
    layers: int = typer.Option(3, "--layers", help="Number of pyramid layers."),
    kappa: float = typer.Option(12.0, "--kappa", help="Inlier detection kappa."),
    num_neighbours: int = typer.Option(3, "--num-neighbours", help="Correspondence neighbours."),
) -> None:
    """Run pyramid (multi-resolution nonrigid) mesh registration on FLOATING → TARGET."""
    trimesh = _require_trimesh()
    import meshmonk  # noqa: PLC0415

    floating_mesh = _load_mesh(trimesh, floating)
    target_mesh = _load_mesh(trimesh, target)

    result = meshmonk.pyramid_register(
        floating=floating_mesh,
        target=target_mesh,
        num_iterations=iterations,
        num_pyramid_layers=layers,
        inlier_kappa=kappa,
        correspondences_num_neighbours=num_neighbours,
    )

    _save_mesh(trimesh, result.aligned_vertices, floating_mesh.faces, out)
    typer.echo(f"Saved result to {out}")


# ---------------------------------------------------------------------------
# demo subcommand
# ---------------------------------------------------------------------------

@app.command()
def demo(
    mode: str = typer.Argument("rigid", help="Registration mode: rigid, nonrigid, or pyramid."),
    download: bool = typer.Option(False, "--download", help="Download demo meshes to ~/.cache/meshmonk/."),
) -> None:
    """Download demo meshes and/or run a demo registration.

    Without --download: looks for meshes in ~/.cache/meshmonk/ then
    falls back to <repo-root>/data/ if running from the repository.

    With --download: fetches Template.obj and demoFace.obj to ~/.cache/meshmonk/.
    """
    if download:
        _demo_download()
        return

    _demo_run(mode)


def _demo_download() -> None:
    """Fetch demo meshes to the cache directory."""
    if _TEMPLATE_URL is None or _DEMO_FACE_URL is None:
        typer.echo("TODO: set download URLs")
        raise typer.Exit(code=0)

    _CACHE_DIR.mkdir(parents=True, exist_ok=True)

    try:
        import urllib.request  # noqa: PLC0415
    except ImportError:
        typer.echo("urllib.request not available — cannot download.", err=True)
        raise typer.Exit(code=1)

    for url, name in [(_TEMPLATE_URL, "Template.obj"), (_DEMO_FACE_URL, "demoFace.obj")]:
        dest = _CACHE_DIR / name
        typer.echo(f"Downloading {name} …")
        urllib.request.urlretrieve(url, dest)  # noqa: S310
        typer.echo(f"  → {dest}")


def _demo_run(mode: str) -> None:
    """Run a demo registration in *mode* using available meshes."""
    paths = _find_demo_meshes()
    if paths is None:
        typer.echo(
            "Demo meshes not found.  Run with --download to fetch them:\n"
            "  meshmonk demo --download\n\n"
            "Or, if you have the repository checked out with data/ present,\n"
            "run from the repository root directory."
        )
        raise typer.Exit(code=1)

    template_path, demo_face_path = paths
    _CACHE_DIR.mkdir(parents=True, exist_ok=True)
    out_path = _CACHE_DIR / "demo_result.obj"

    trimesh = _require_trimesh()
    import meshmonk  # noqa: PLC0415

    floating_mesh = _load_mesh(trimesh, template_path)
    target_mesh = _load_mesh(trimesh, demo_face_path)

    if mode == "rigid":
        result = meshmonk.rigid_register(floating=floating_mesh, target=target_mesh)
    elif mode == "nonrigid":
        result = meshmonk.nonrigid_register(floating=floating_mesh, target=target_mesh)
    elif mode == "pyramid":
        result = meshmonk.pyramid_register(floating=floating_mesh, target=target_mesh)
    else:
        typer.echo(f"Unknown mode {mode!r}.  Choose: rigid, nonrigid, pyramid.", err=True)
        raise typer.Exit(code=1)

    _save_mesh(trimesh, result.aligned_vertices, floating_mesh.faces, str(out_path))
    typer.echo(f"Demo result saved to {out_path}")
