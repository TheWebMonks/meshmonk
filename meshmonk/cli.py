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

DEMO_ASSETS = {
    "Template.obj": {
        "url": "https://github.com/jsnyde0/meshmonk/releases/latest/download/Template.obj",
        "sha256": None,  # populated at first release — set to None until v0.3.0 is tagged
    },
    "demoFace.obj": {
        "url": "https://github.com/jsnyde0/meshmonk/releases/latest/download/demoFace.obj",
        "sha256": None,  # redistribution rights TBD — see Prerequisites in design doc
    },
}

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
    import hashlib
    import urllib.error
    import urllib.request

    _CACHE_DIR.mkdir(parents=True, exist_ok=True)

    for name, asset in DEMO_ASSETS.items():
        dest = _CACHE_DIR / name
        typer.echo(f"Downloading {name} ...")
        try:
            with urllib.request.urlopen(asset["url"], timeout=30) as resp:  # noqa: S310
                dest.write_bytes(resp.read())
        except (urllib.error.URLError, OSError) as e:
            dest.unlink(missing_ok=True)
            typer.echo(str(e), err=True)
            raise typer.Exit(code=1)
        if asset["sha256"] is not None:
            digest = hashlib.sha256(dest.read_bytes()).hexdigest()
            if digest != asset["sha256"]:
                dest.unlink(missing_ok=True)
                typer.echo(
                    f"SHA-256 mismatch for {name}. "
                    "Your meshmonk version may be out of date. "
                    "Run: pip install --upgrade meshmonk",
                    err=True,
                )
                raise typer.Exit(code=1)
        else:
            typer.echo(f"  Warning: integrity check not configured for {name}")
        typer.echo(f"  -> {dest}")


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
