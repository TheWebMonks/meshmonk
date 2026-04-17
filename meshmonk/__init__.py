"""MeshMonk — Python-first 3D mesh registration library."""
from __future__ import annotations

__version__ = "0.1.0"

import warnings
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

# ---------------------------------------------------------------------------
# C++ extension module
# ---------------------------------------------------------------------------
from meshmonk._meshmonk_core import (  # noqa: E402
    # Enums
    LogLevel,
    RegistrationError,
    # Param structs
    CorrespondenceParams,
    DownsampleSchedule,
    InlierParams,
    NonrigidParams,
    PyramidParams,
    RigidParams,
    ViscoElasticParams,
    # Transform
    RigidTransform,
    # Result structs (raw C++ objects, wrapped by Python dataclasses below)
    RigidResult as _RigidResult,
    NonrigidResult as _NonrigidResult,
    PyramidResult as _PyramidResult,
    # High-level pipelines (raw)
    rigid_registration as _rigid_registration,
    nonrigid_registration as _nonrigid_registration,
    pyramid_registration as _pyramid_registration,
    # Low-level primitives — re-exported directly
    compute_correspondences,
    compute_inlier_weights,
    compute_rigid_transform,
    compute_nonrigid_transform,
    downsample_mesh,
    scale_shift_mesh,
    compute_normals,
    set_log_level as _set_log_level_core,
)

if TYPE_CHECKING:
    from pathlib import Path

# ---------------------------------------------------------------------------
# Public re-exports for primitives
# ---------------------------------------------------------------------------
__all__ = [
    "__version__",
    # Enums
    "LogLevel",
    "RegistrationError",
    # Param structs
    "CorrespondenceParams",
    "DownsampleSchedule",
    "InlierParams",
    "NonrigidParams",
    "PyramidParams",
    "RigidParams",
    "ViscoElasticParams",
    # Transform
    "RigidTransform",
    # Python result dataclasses
    "RigidRegResult",
    "NonrigidRegResult",
    "PyramidRegResult",
    # High-level registration functions
    "rigid_register",
    "nonrigid_register",
    "pyramid_register",
    # Helpers
    "features_from_vertices",
    "set_log_level",
    # Low-level primitives
    "compute_correspondences",
    "compute_inlier_weights",
    "compute_rigid_transform",
    "compute_nonrigid_transform",
    "downsample_mesh",
    "scale_shift_mesh",
    "compute_normals",
]


# ---------------------------------------------------------------------------
# Python result dataclasses
# ---------------------------------------------------------------------------

@dataclass
class RigidRegResult:
    """Result of rigid_register()."""
    aligned_features: np.ndarray   # (N, 6) float32 — positions + normals
    transform: RigidTransform       # 4x4 SE(3) rigid transform
    iterations_run: int

    @property
    def aligned_vertices(self) -> np.ndarray:
        """(N, 3) float32 — aligned positions only."""
        return self.aligned_features[:, :3]


@dataclass
class NonrigidRegResult:
    """Result of nonrigid_register()."""
    aligned_features: np.ndarray      # (N, 6) float32
    final_inlier_weights: np.ndarray  # (N,) float32
    displacement_field: np.ndarray    # (N, 3) float32
    iterations_run: int

    @property
    def aligned_vertices(self) -> np.ndarray:
        """(N, 3) float32 — aligned positions only."""
        return self.aligned_features[:, :3]


@dataclass
class PyramidRegResult:
    """Result of pyramid_register()."""
    aligned_features: np.ndarray      # (N, 6) float32
    final_inlier_weights: np.ndarray  # (N,) float32
    displacement_field: np.ndarray    # (N, 3) float32
    per_layer_iterations: list[int]

    @property
    def aligned_vertices(self) -> np.ndarray:
        """(N, 3) float32 — aligned positions only."""
        return self.aligned_features[:, :3]


# ---------------------------------------------------------------------------
# set_log_level convenience wrapper
# ---------------------------------------------------------------------------

_LOG_LEVEL_MAP = {
    "silent":  LogLevel.Silent,
    "error":   LogLevel.Error,
    "warning": LogLevel.Warning,
    "info":    LogLevel.Info,
    "debug":   LogLevel.Debug,
}


def set_log_level(level: str) -> None:
    """Set the C++ log level.

    Parameters
    ----------
    level:
        One of 'silent', 'error', 'warning', 'info', 'debug' (case-insensitive).
    """
    key = level.lower()
    if key not in _LOG_LEVEL_MAP:
        raise ValueError(
            f"Unknown log level {level!r}. "
            f"Valid values: {sorted(_LOG_LEVEL_MAP.keys())}"
        )
    _set_log_level_core(_LOG_LEVEL_MAP[key])


# ---------------------------------------------------------------------------
# features_from_vertices helper
# ---------------------------------------------------------------------------

def features_from_vertices(V: np.ndarray, F: np.ndarray) -> np.ndarray:
    """Build a (N, 6) feature matrix from positions and face topology.

    Columns 0-2 are float32 positions from *V*; columns 3-5 are face-area-
    weighted vertex normals computed via compute_normals(V, F).

    Parameters
    ----------
    V:
        (N, 3) float32 array of vertex positions.
    F:
        (M, 3) int32 array of triangular faces.

    Returns
    -------
    np.ndarray
        (N, 6) float32 feature matrix.
    """
    V = np.asarray(V, dtype=np.float32)
    F = np.asarray(F, dtype=np.int32)
    normals = np.asarray(compute_normals(V, F), dtype=np.float32)
    features = np.empty((V.shape[0], 6), dtype=np.float32)
    features[:, :3] = V
    features[:, 3:] = normals
    return features


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _load_trimesh(path):
    """Load a mesh using trimesh (requires meshmonk[io])."""
    try:
        import trimesh  # noqa: PLC0415
    except ImportError as exc:
        raise ImportError(
            "Loading meshes from file paths requires the 'io' extra. "
            "Install with: pip install 'meshmonk[io]'"
        ) from exc
    return trimesh.load(str(path))


def _mesh_to_arrays(mesh_or_path, normals_override=None, force_recompute=False):
    """Convert a duck-typed mesh object or path to (features, faces, flags).

    Parameters
    ----------
    mesh_or_path:
        trimesh.Trimesh, str, or pathlib.Path.
    normals_override:
        Optional explicit (N, 3) normals array.
    force_recompute:
        If True, recompute normals even when mesh.vertex_normals is available.

    Returns
    -------
    (features, faces, flags)
        features: (N, 6) float32
        faces:    (M, 3) int32
        flags:    (N,)   float32 (all ones)
    """
    import pathlib  # noqa: PLC0415

    if isinstance(mesh_or_path, (str, pathlib.Path)):
        mesh = _load_trimesh(mesh_or_path)
    else:
        mesh = mesh_or_path

    V = np.asarray(mesh.vertices, dtype=np.float32)
    # trimesh returns int64 faces — silently cast to int32 at Python boundary
    F = np.asarray(mesh.faces, dtype=np.int32)
    flags = np.ones(V.shape[0], dtype=np.float32)

    if normals_override is not None:
        N = np.asarray(normals_override, dtype=np.float32)
    elif (not force_recompute
          and hasattr(mesh, "vertex_normals")
          and mesh.vertex_normals is not None):
        vn = np.asarray(mesh.vertex_normals, dtype=np.float32)
        norms = np.linalg.norm(vn, axis=1)
        if np.any(norms > 0):
            N = vn
        else:
            N = np.asarray(compute_normals(V, F), dtype=np.float32)
    else:
        N = np.asarray(compute_normals(V, F), dtype=np.float32)

    features = np.empty((V.shape[0], 6), dtype=np.float32)
    features[:, :3] = V
    features[:, 3:] = N
    return features, F, flags


def _check_normals_and_warn(features: np.ndarray) -> np.ndarray:
    """If normals (cols 3-5) are all-zero, recompute from positions and warn.

    This requires faces to be available, but for the raw-array path we don't
    have faces here. The check and recompute happens at the caller level.
    """
    return features


def _apply_rigid_kwargs(params: RigidParams, kwargs: dict) -> RigidParams:
    """Apply kwarg overrides to a RigidParams struct."""
    if "correspondences_symmetric" in kwargs:
        params.correspondences.symmetric = kwargs["correspondences_symmetric"]
    if "correspondences_num_neighbours" in kwargs:
        params.correspondences.num_neighbours = kwargs["correspondences_num_neighbours"]
    if "correspondences_flag_threshold" in kwargs:
        params.correspondences.flag_threshold = kwargs["correspondences_flag_threshold"]
    if "correspondences_equalize_push_pull" in kwargs:
        params.correspondences.equalize_push_pull = kwargs["correspondences_equalize_push_pull"]
    if "inlier_kappa" in kwargs:
        params.inliers.kappa = kwargs["inlier_kappa"]
    if "inlier_use_orientation" in kwargs:
        params.inliers.use_orientation = kwargs["inlier_use_orientation"]
    if "num_iterations" in kwargs:
        params.num_iterations = kwargs["num_iterations"]
    if "use_scaling" in kwargs:
        params.use_scaling = kwargs["use_scaling"]
    return params


def _apply_nonrigid_kwargs(params: NonrigidParams, kwargs: dict) -> NonrigidParams:
    """Apply kwarg overrides to a NonrigidParams struct."""
    if "correspondences_symmetric" in kwargs:
        params.correspondences.symmetric = kwargs["correspondences_symmetric"]
    if "correspondences_num_neighbours" in kwargs:
        params.correspondences.num_neighbours = kwargs["correspondences_num_neighbours"]
    if "correspondences_flag_threshold" in kwargs:
        params.correspondences.flag_threshold = kwargs["correspondences_flag_threshold"]
    if "correspondences_equalize_push_pull" in kwargs:
        params.correspondences.equalize_push_pull = kwargs["correspondences_equalize_push_pull"]
    if "inlier_kappa" in kwargs:
        params.inliers.kappa = kwargs["inlier_kappa"]
    if "inlier_use_orientation" in kwargs:
        params.inliers.use_orientation = kwargs["inlier_use_orientation"]
    if "transform_sigma" in kwargs:
        params.transform.sigma = kwargs["transform_sigma"]
    if "transform_num_viscous_iterations_start" in kwargs:
        params.transform.num_viscous_iterations_start = kwargs["transform_num_viscous_iterations_start"]
    if "transform_num_viscous_iterations_end" in kwargs:
        params.transform.num_viscous_iterations_end = kwargs["transform_num_viscous_iterations_end"]
    if "transform_num_elastic_iterations_start" in kwargs:
        params.transform.num_elastic_iterations_start = kwargs["transform_num_elastic_iterations_start"]
    if "transform_num_elastic_iterations_end" in kwargs:
        params.transform.num_elastic_iterations_end = kwargs["transform_num_elastic_iterations_end"]
    if "num_iterations" in kwargs:
        params.num_iterations = kwargs["num_iterations"]
    return params


def _apply_pyramid_kwargs(params: PyramidParams, kwargs: dict, explicit_kwargs: set) -> PyramidParams:
    """Apply kwarg overrides to a PyramidParams struct.

    Also applies the MATLAB convention: viscous/elastic start = num_iterations
    when the user has not explicitly set them.
    """
    if "correspondences_symmetric" in kwargs:
        params.correspondences.symmetric = kwargs["correspondences_symmetric"]
    if "correspondences_num_neighbours" in kwargs:
        params.correspondences.num_neighbours = kwargs["correspondences_num_neighbours"]
    if "correspondences_flag_threshold" in kwargs:
        params.correspondences.flag_threshold = kwargs["correspondences_flag_threshold"]
    if "correspondences_equalize_push_pull" in kwargs:
        params.correspondences.equalize_push_pull = kwargs["correspondences_equalize_push_pull"]
    if "inlier_kappa" in kwargs:
        params.inliers.kappa = kwargs["inlier_kappa"]
    if "inlier_use_orientation" in kwargs:
        params.inliers.use_orientation = kwargs["inlier_use_orientation"]
    if "transform_sigma" in kwargs:
        params.transform.sigma = kwargs["transform_sigma"]
    if "transform_num_viscous_iterations_start" in kwargs:
        params.transform.num_viscous_iterations_start = kwargs["transform_num_viscous_iterations_start"]
    if "transform_num_viscous_iterations_end" in kwargs:
        params.transform.num_viscous_iterations_end = kwargs["transform_num_viscous_iterations_end"]
    if "transform_num_elastic_iterations_start" in kwargs:
        params.transform.num_elastic_iterations_start = kwargs["transform_num_elastic_iterations_start"]
    if "transform_num_elastic_iterations_end" in kwargs:
        params.transform.num_elastic_iterations_end = kwargs["transform_num_elastic_iterations_end"]
    if "downsample_float_start" in kwargs:
        params.downsample.float_start = kwargs["downsample_float_start"]
    if "downsample_target_start" in kwargs:
        params.downsample.target_start = kwargs["downsample_target_start"]
    if "downsample_float_end" in kwargs:
        params.downsample.float_end = kwargs["downsample_float_end"]
    if "downsample_target_end" in kwargs:
        params.downsample.target_end = kwargs["downsample_target_end"]
    if "num_iterations" in kwargs:
        params.num_iterations = kwargs["num_iterations"]
    if "num_pyramid_layers" in kwargs:
        params.num_pyramid_layers = kwargs["num_pyramid_layers"]

    # MATLAB convention: auto-populate viscous/elastic start = num_iterations
    # when the user has not explicitly set them
    if "transform_num_viscous_iterations_start" not in explicit_kwargs:
        params.transform.num_viscous_iterations_start = params.num_iterations
    if "transform_num_elastic_iterations_start" not in explicit_kwargs:
        params.transform.num_elastic_iterations_start = params.num_iterations

    return params


# ---------------------------------------------------------------------------
# High-level registration functions
# ---------------------------------------------------------------------------

def rigid_register(
    floating=None,
    target=None,
    *,
    floating_features=None,
    target_features=None,
    floating_faces=None,
    target_faces=None,
    floating_flags=None,
    target_flags=None,
    normals=None,
    compute_normals_flag=False,
    **kwargs,
) -> RigidRegResult:
    """Run rigid (SE(3)) mesh registration.

    Two call patterns are supported:

    Pattern A — duck-typed mesh objects or path strings::

        result = meshmonk.rigid_register(floating=mesh_obj, target="demoFace.obj")

    Pattern B — raw numpy arrays::

        result = meshmonk.rigid_register(
            floating_features=V_float,
            target_features=V_target,
            floating_faces=F_float,
            target_faces=F_target,
            floating_flags=flags_float,
            target_flags=flags_target,
        )

    Parameters
    ----------
    floating, target:
        trimesh.Trimesh objects, file path strings, or pathlib.Path objects.
    floating_features, target_features:
        (N, 6) float32 feature matrices (positions + normals).
    floating_faces, target_faces:
        (M, 3) int32 face index arrays.
    floating_flags, target_flags:
        (N,) float32 per-vertex flags (1.0 = active, 0.0 = masked).
    normals:
        Optional explicit (N, 3) normals for Pattern A.
    compute_normals_flag:
        If True, force recomputation of normals even when mesh.vertex_normals is available.
    **kwargs:
        Param overrides. See kwarg flattening rule in design doc.

    Returns
    -------
    RigidRegResult
    """
    if floating is not None or target is not None:
        # Pattern A: duck-typed mesh objects
        feat_float, faces_float, flags_float = _mesh_to_arrays(
            floating, normals_override=normals, force_recompute=compute_normals_flag
        )
        feat_target, faces_target, flags_target = _mesh_to_arrays(
            target, force_recompute=compute_normals_flag
        )
    else:
        # Pattern B: raw numpy arrays
        feat_float  = np.asarray(floating_features, dtype=np.float32)
        feat_target = np.asarray(target_features, dtype=np.float32)
        faces_float  = np.asarray(floating_faces, dtype=np.int32)
        faces_target = np.asarray(target_faces, dtype=np.int32)
        flags_float  = np.asarray(floating_flags, dtype=np.float32)
        flags_target = np.asarray(target_flags, dtype=np.float32)

        # Auto-recompute normals if all-zero
        if np.all(feat_float[:, 3:] == 0.0):
            warnings.warn(
                "floating_features normals are all-zero; recomputing.",
                stacklevel=2,
            )
            normals_recomputed = np.asarray(
                compute_normals(feat_float[:, :3], faces_float), dtype=np.float32
            )
            feat_float = feat_float.copy()
            feat_float[:, 3:] = normals_recomputed

    params = RigidParams()
    params = _apply_rigid_kwargs(params, kwargs)

    raw = _rigid_registration(
        feat_float, feat_target,
        faces_float, faces_target,
        flags_float, flags_target,
        params,
    )
    return RigidRegResult(
        aligned_features=np.asarray(raw.aligned_features, dtype=np.float32),
        transform=raw.transform,
        iterations_run=raw.iterations_run,
    )


def nonrigid_register(
    floating=None,
    target=None,
    *,
    floating_features=None,
    target_features=None,
    floating_faces=None,
    target_faces=None,
    floating_flags=None,
    target_flags=None,
    normals=None,
    compute_normals_flag=False,
    **kwargs,
) -> NonrigidRegResult:
    """Run nonrigid (viscoelastic) mesh registration.

    Accepts the same two call patterns as rigid_register().

    Returns
    -------
    NonrigidRegResult
    """
    if floating is not None or target is not None:
        feat_float, faces_float, flags_float = _mesh_to_arrays(
            floating, normals_override=normals, force_recompute=compute_normals_flag
        )
        feat_target, faces_target, flags_target = _mesh_to_arrays(
            target, force_recompute=compute_normals_flag
        )
    else:
        feat_float  = np.asarray(floating_features, dtype=np.float32)
        feat_target = np.asarray(target_features, dtype=np.float32)
        faces_float  = np.asarray(floating_faces, dtype=np.int32)
        faces_target = np.asarray(target_faces, dtype=np.int32)
        flags_float  = np.asarray(floating_flags, dtype=np.float32)
        flags_target = np.asarray(target_flags, dtype=np.float32)

        if np.all(feat_float[:, 3:] == 0.0):
            warnings.warn(
                "floating_features normals are all-zero; recomputing.",
                stacklevel=2,
            )
            normals_recomputed = np.asarray(
                compute_normals(feat_float[:, :3], faces_float), dtype=np.float32
            )
            feat_float = feat_float.copy()
            feat_float[:, 3:] = normals_recomputed

    params = NonrigidParams()
    params = _apply_nonrigid_kwargs(params, kwargs)

    raw = _nonrigid_registration(
        feat_float, feat_target,
        faces_float, faces_target,
        flags_float, flags_target,
        params,
    )
    return NonrigidRegResult(
        aligned_features=np.asarray(raw.aligned_features, dtype=np.float32),
        final_inlier_weights=np.asarray(raw.final_inlier_weights, dtype=np.float32),
        displacement_field=np.asarray(raw.displacement_field, dtype=np.float32),
        iterations_run=raw.iterations_run,
    )


def pyramid_register(
    floating=None,
    target=None,
    *,
    floating_features=None,
    target_features=None,
    floating_faces=None,
    target_faces=None,
    floating_flags=None,
    target_flags=None,
    normals=None,
    compute_normals_flag=False,
    **kwargs,
) -> PyramidRegResult:
    """Run pyramid (multi-resolution) nonrigid mesh registration.

    Accepts the same two call patterns as rigid_register().

    Special behaviour: viscous/elastic start iterations are auto-set to
    num_iterations when not explicitly passed (matching MATLAB convention).

    Returns
    -------
    PyramidRegResult
    """
    if floating is not None or target is not None:
        feat_float, faces_float, flags_float = _mesh_to_arrays(
            floating, normals_override=normals, force_recompute=compute_normals_flag
        )
        feat_target, faces_target, flags_target = _mesh_to_arrays(
            target, force_recompute=compute_normals_flag
        )
    else:
        feat_float  = np.asarray(floating_features, dtype=np.float32)
        feat_target = np.asarray(target_features, dtype=np.float32)
        faces_float  = np.asarray(floating_faces, dtype=np.int32)
        faces_target = np.asarray(target_faces, dtype=np.int32)
        flags_float  = np.asarray(floating_flags, dtype=np.float32)
        flags_target = np.asarray(target_flags, dtype=np.float32)

        if np.all(feat_float[:, 3:] == 0.0):
            warnings.warn(
                "floating_features normals are all-zero; recomputing.",
                stacklevel=2,
            )
            normals_recomputed = np.asarray(
                compute_normals(feat_float[:, :3], faces_float), dtype=np.float32
            )
            feat_float = feat_float.copy()
            feat_float[:, 3:] = normals_recomputed

    # Pass the set of explicitly provided kwargs so _apply_pyramid_kwargs
    # knows which viscous/elastic start values to auto-populate
    params = PyramidParams()
    # Override flag_threshold to MATLAB pyramid default (0.999) — Python
    # wrapper responsibility, not the C++ struct default
    params.correspondences.flag_threshold = 0.999
    params = _apply_pyramid_kwargs(params, kwargs, explicit_kwargs=set(kwargs.keys()))

    raw = _pyramid_registration(
        feat_float, feat_target,
        faces_float, faces_target,
        flags_float, flags_target,
        params,
    )
    return PyramidRegResult(
        aligned_features=np.asarray(raw.aligned_features, dtype=np.float32),
        final_inlier_weights=np.asarray(raw.final_inlier_weights, dtype=np.float32),
        displacement_field=np.asarray(raw.displacement_field, dtype=np.float32),
        per_layer_iterations=list(raw.per_layer_iterations),
    )
