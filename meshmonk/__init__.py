"""MeshMonk — Python-first 3D mesh registration library."""

from __future__ import annotations

from typing import Any, Mapping
from typing_extensions import Unpack
from meshmonk._types import RigidKwargs, NonrigidKwargs, PyramidKwargs

try:
    from importlib.metadata import version

    __version__ = version("meshmonk")
except Exception:
    __version__ = "0.0.0.dev0"  # fallback for uninstalled dev runs

from dataclasses import dataclass

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
    # Exception
    MeshMonkError,
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

# ---------------------------------------------------------------------------
# Add .code property to MeshMonkError for programmatic error dispatch
# ---------------------------------------------------------------------------


def _meshmonk_error_code(self):
    """The RegistrationError code for this error, or None."""
    raw = getattr(self, "_code", None)
    if raw is not None:
        return RegistrationError(raw)
    return None


MeshMonkError.code = property(_meshmonk_error_code)  # pyright: ignore[reportAttributeAccessIssue]

# ---------------------------------------------------------------------------
# Public re-exports for primitives
# ---------------------------------------------------------------------------
__all__ = [
    "__version__",
    # Enums
    "LogLevel",
    "RegistrationError",
    # Exception
    "MeshMonkError",
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

    aligned_features: np.ndarray  # (N, 6) float32 — positions + normals
    transform: RigidTransform  # 4x4 SE(3) rigid transform
    iterations_run: int

    @property
    def aligned_vertices(self) -> np.ndarray:
        """(N, 3) float32 — aligned positions only."""
        return self.aligned_features[:, :3]


@dataclass
class NonrigidRegResult:
    """Result of nonrigid_register()."""

    aligned_features: np.ndarray  # (N, 6) float32
    final_inlier_weights: np.ndarray  # (N,) float32
    displacement_field: np.ndarray  # (N, 3) float32
    iterations_run: int

    @property
    def aligned_vertices(self) -> np.ndarray:
        """(N, 3) float32 — aligned positions only."""
        return self.aligned_features[:, :3]


@dataclass
class PyramidRegResult:
    """Result of pyramid_register()."""

    aligned_features: np.ndarray  # (N, 6) float32
    final_inlier_weights: np.ndarray  # (N,) float32
    displacement_field: np.ndarray  # (N, 3) float32
    per_layer_iterations: list[int]

    @property
    def aligned_vertices(self) -> np.ndarray:
        """(N, 3) float32 — aligned positions only."""
        return self.aligned_features[:, :3]


# ---------------------------------------------------------------------------
# set_log_level convenience wrapper
# ---------------------------------------------------------------------------

_LOG_LEVEL_MAP = {
    "silent": LogLevel.Silent,
    "error": LogLevel.Error,
    "warning": LogLevel.Warning,
    "info": LogLevel.Info,
    "debug": LogLevel.Debug,
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
    from typing import cast  # noqa: PLC0415

    return cast("trimesh.Trimesh", trimesh.load(str(path)))


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
    _user_flags = getattr(mesh, "flags", None)
    if _user_flags is not None:
        flags = np.asarray(_user_flags, dtype=np.float32)
        if flags.shape != (V.shape[0],):
            raise ValueError(
                f"mesh.flags shape {flags.shape} does not match vertex count ({V.shape[0]},)"
            )

    if normals_override is not None:
        N = np.asarray(normals_override, dtype=np.float32)
        if N.shape != (V.shape[0], 3):
            raise ValueError(
                f"normals_override shape {N.shape} does not match expected ({V.shape[0]}, 3)"
            )
    elif (
        not force_recompute
        and hasattr(mesh, "vertex_normals")
        and mesh.vertex_normals is not None
    ):
        vn = np.asarray(mesh.vertex_normals, dtype=np.float32)
        norms = np.linalg.norm(vn, axis=1)
        if np.all(norms > 0):
            N = vn
        else:
            N = np.asarray(compute_normals(V, F), dtype=np.float32)
    else:
        N = np.asarray(compute_normals(V, F), dtype=np.float32)

    features = np.empty((V.shape[0], 6), dtype=np.float32)
    features[:, :3] = V
    features[:, 3:] = N
    return features, F, flags


_RIGID_KNOWN_KWARGS = {
    "correspondences_symmetric",
    "correspondences_num_neighbours",
    "correspondences_flag_threshold",
    "correspondences_equalize_push_pull",
    "inlier_kappa",
    "inlier_use_orientation",
    "num_iterations",
    "use_scaling",
}

_NONRIGID_KNOWN_KWARGS = {
    "correspondences_symmetric",
    "correspondences_num_neighbours",
    "correspondences_flag_threshold",
    "correspondences_equalize_push_pull",
    "inlier_kappa",
    "inlier_use_orientation",
    "transform_sigma",
    "transform_num_viscous_iterations_start",
    "transform_num_viscous_iterations_end",
    "transform_num_elastic_iterations_start",
    "transform_num_elastic_iterations_end",
    "num_iterations",
}

_PYRAMID_KNOWN_KWARGS = {
    "correspondences_symmetric",
    "correspondences_num_neighbours",
    "correspondences_flag_threshold",
    "correspondences_equalize_push_pull",
    "inlier_kappa",
    "inlier_use_orientation",
    "transform_sigma",
    "transform_num_viscous_iterations_start",
    "transform_num_viscous_iterations_end",
    "transform_num_elastic_iterations_start",
    "transform_num_elastic_iterations_end",
    "downsample_float_start",
    "downsample_target_start",
    "downsample_float_end",
    "downsample_target_end",
    "num_iterations",
    "num_pyramid_layers",
}


def _apply_shared_kwargs(params, kwargs: Mapping[str, Any]) -> None:
    """Apply correspondence + inlier kwarg overrides shared by all registration types."""
    if "correspondences_symmetric" in kwargs:
        params.correspondences.symmetric = kwargs["correspondences_symmetric"]
    if "correspondences_num_neighbours" in kwargs:
        params.correspondences.num_neighbours = kwargs["correspondences_num_neighbours"]
    if "correspondences_flag_threshold" in kwargs:
        params.correspondences.flag_threshold = kwargs["correspondences_flag_threshold"]
    if "correspondences_equalize_push_pull" in kwargs:
        params.correspondences.equalize_push_pull = kwargs[
            "correspondences_equalize_push_pull"
        ]
    if "inlier_kappa" in kwargs:
        params.inliers.kappa = kwargs["inlier_kappa"]
    if "inlier_use_orientation" in kwargs:
        params.inliers.use_orientation = kwargs["inlier_use_orientation"]


def _apply_transform_kwargs(params, kwargs: Mapping[str, Any]) -> None:
    """Apply transform kwarg overrides shared by nonrigid and pyramid registration."""
    if "transform_sigma" in kwargs:
        params.transform.sigma = kwargs["transform_sigma"]
    if "transform_num_viscous_iterations_start" in kwargs:
        params.transform.num_viscous_iterations_start = kwargs[
            "transform_num_viscous_iterations_start"
        ]
    if "transform_num_viscous_iterations_end" in kwargs:
        params.transform.num_viscous_iterations_end = kwargs[
            "transform_num_viscous_iterations_end"
        ]
    if "transform_num_elastic_iterations_start" in kwargs:
        params.transform.num_elastic_iterations_start = kwargs[
            "transform_num_elastic_iterations_start"
        ]
    if "transform_num_elastic_iterations_end" in kwargs:
        params.transform.num_elastic_iterations_end = kwargs[
            "transform_num_elastic_iterations_end"
        ]


def _apply_rigid_kwargs(params: RigidParams, kwargs: Mapping[str, Any]) -> RigidParams:
    """Apply kwarg overrides to a RigidParams struct."""
    unknown = set(kwargs) - _RIGID_KNOWN_KWARGS
    if unknown:
        raise TypeError(
            f"rigid_register() got unexpected keyword arguments: {sorted(unknown)}"
        )
    _apply_shared_kwargs(params, kwargs)
    if "num_iterations" in kwargs:
        params.num_iterations = kwargs["num_iterations"]
    if "use_scaling" in kwargs:
        params.use_scaling = kwargs["use_scaling"]
    return params


def _apply_nonrigid_kwargs(params: NonrigidParams, kwargs: Mapping[str, Any]) -> NonrigidParams:
    """Apply kwarg overrides to a NonrigidParams struct."""
    unknown = set(kwargs) - _NONRIGID_KNOWN_KWARGS
    if unknown:
        raise TypeError(
            f"nonrigid_register() got unexpected keyword arguments: {sorted(unknown)}"
        )
    _apply_shared_kwargs(params, kwargs)
    _apply_transform_kwargs(params, kwargs)
    if "num_iterations" in kwargs:
        params.num_iterations = kwargs["num_iterations"]
    return params


def _apply_pyramid_kwargs(
    params: PyramidParams, kwargs: Mapping[str, Any], explicit_kwargs: set[str]
) -> PyramidParams:
    """Apply kwarg overrides to a PyramidParams struct.

    Also applies the MATLAB convention: viscous/elastic start = num_iterations
    when the user has not explicitly set them.
    """
    unknown = set(kwargs) - _PYRAMID_KNOWN_KWARGS
    if unknown:
        raise TypeError(
            f"pyramid_register() got unexpected keyword arguments: {sorted(unknown)}"
        )
    _apply_shared_kwargs(params, kwargs)
    _apply_transform_kwargs(params, kwargs)
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
# Shared array-preparation helper
# ---------------------------------------------------------------------------


def _prepare_arrays(
    floating=None,
    target=None,
    *,
    floating_features=None,
    target_features=None,
    floating_faces=None,
    target_faces=None,
    floating_normals_override=None,
    force_recompute_normals=False,
    floating_flags=None,
    target_flags=None,
):
    """Validate inputs and convert to (feat_float, feat_target, faces_float, faces_target, flags_float, flags_target).

    Supports both Pattern A (duck-typed mesh objects / file paths) and
    Pattern B (raw numpy arrays).  In Pattern B, flags default to all-ones
    when not provided.
    """
    import warnings as _warnings
    from meshmonk._meshmonk_core import compute_normals as _compute_normals

    # Symmetry check
    if (floating is None) != (target is None):
        raise ValueError("both floating and target must be provided, or neither")
    # Mixed Pattern A + B check
    pattern_b_args = any(
        x is not None
        for x in [
            floating_features,
            target_features,
            floating_faces,
            target_faces,
            floating_flags,
            target_flags,
        ]
    )
    if floating is not None and pattern_b_args:
        raise TypeError(
            "Cannot mix Pattern A (floating/target) with Pattern B "
            "(floating_features/target_features) arguments"
        )
    pattern_a = floating is not None
    if pattern_a:
        feat_float, faces_float, flags_float = _mesh_to_arrays(
            floating,
            normals_override=floating_normals_override,
            force_recompute=force_recompute_normals,
        )
        feat_target, faces_target, flags_target = _mesh_to_arrays(
            target, force_recompute=force_recompute_normals
        )
        if floating_flags is not None:
            flags_float = np.asarray(floating_flags, dtype=np.float32)
        if target_flags is not None:
            flags_target = np.asarray(target_flags, dtype=np.float32)
    else:
        required = {
            "floating_features": floating_features,
            "target_features": target_features,
            "floating_faces": floating_faces,
            "target_faces": target_faces,
        }
        missing = [k for k, v in required.items() if v is None]
        if missing:
            raise ValueError(
                "Pattern B requires all array arguments. Missing: " + ", ".join(missing)
            )
        feat_float = np.asarray(floating_features, dtype=np.float32)
        feat_target = np.asarray(target_features, dtype=np.float32)
        faces_float = np.asarray(floating_faces, dtype=np.int32)
        faces_target = np.asarray(target_faces, dtype=np.int32)
        n_float = feat_float.shape[0]
        n_target = feat_target.shape[0]
        flags_float = (
            np.asarray(floating_flags, dtype=np.float32)
            if floating_flags is not None
            else np.ones(n_float, dtype=np.float32)
        )
        flags_target = (
            np.asarray(target_flags, dtype=np.float32)
            if target_flags is not None
            else np.ones(n_target, dtype=np.float32)
        )
        # Auto-recompute normals if all-zero
        if feat_float.shape[1] == 6 and np.all(feat_float[:, 3:] == 0.0):
            _warnings.warn(
                "floating_features normals are all-zero; recomputing.", stacklevel=3
            )
            normals_r = np.asarray(
                _compute_normals(feat_float[:, :3], faces_float), dtype=np.float32
            )
            feat_float = feat_float.copy()
            feat_float[:, 3:] = normals_r
        if feat_target.shape[1] == 6 and np.all(feat_target[:, 3:] == 0.0):
            _warnings.warn(
                "target_features normals are all-zero; recomputing.", stacklevel=3
            )
            normals_r = np.asarray(
                _compute_normals(feat_target[:, :3], faces_target), dtype=np.float32
            )
            feat_target = feat_target.copy()
            feat_target[:, 3:] = normals_r

    return feat_float, feat_target, faces_float, faces_target, flags_float, flags_target


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
    **kwargs: Unpack[RigidKwargs],
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

    Named Parameters
    ----------------
    Pattern A (mesh objects / file paths):
        floating, target : Mesh | None
            trimesh.Trimesh objects, file path strings, or pathlib.Path objects.

    Pattern B (raw numpy arrays):
        floating_features, target_features : NDArray | None
            (N, 6) float32 feature matrices (positions + normals).
        floating_faces, target_faces : NDArray | None
            (M, 3) int32 face index arrays.
        floating_flags, target_flags : NDArray | None
            (N,) float32 per-vertex flags (1.0 = active, 0.0 = masked).

    Pattern B only:
        normals : NDArray | None
            Explicit (N, 3) normals for the floating mesh.
        compute_normals_flag : bool
            Force recomputation of normals even when mesh normals are available.

    nonrigid and pyramid only:
        rigid_params : dict | None
            kwargs for rigid pre-alignment; None skips it.

    **kwargs — per-function scope
    ------------------------------
    The table below lists every accepted kwarg and which functions accept it
    (Y = accepted, — = raises TypeError).

    ==========================================  ================================================  =======  =========  =======
    kwarg                                       params field                                      rigid    nonrigid   pyramid
    ==========================================  ================================================  =======  =========  =======
    correspondences_symmetric                   params.correspondences.symmetric                    Y          Y         Y
    correspondences_num_neighbours              params.correspondences.num_neighbours               Y          Y         Y
    correspondences_flag_threshold [1]          params.correspondences.flag_threshold               Y          Y         Y
    correspondences_equalize_push_pull          params.correspondences.equalize_push_pull           Y          Y         Y
    inlier_kappa                                params.inliers.kappa                                Y          Y         Y
    inlier_use_orientation                      params.inliers.use_orientation                      Y          Y         Y
    num_iterations                              params.num_iterations                               Y          Y         Y
    use_scaling                                 params.use_scaling                                  Y          —         —
    transform_sigma                             params.transform.sigma                              —          Y         Y
    transform_num_viscous_iterations_start [2]  params.transform.num_viscous_iterations_start       —          Y         Y
    transform_num_viscous_iterations_end        params.transform.num_viscous_iterations_end         —          Y         Y
    transform_num_elastic_iterations_start [2]  params.transform.num_elastic_iterations_start       —          Y         Y
    transform_num_elastic_iterations_end        params.transform.num_elastic_iterations_end         —          Y         Y
    downsample_float_start                      params.downsample.float_start                       —          —         Y
    downsample_target_start                     params.downsample.target_start                      —          —         Y
    downsample_float_end                        params.downsample.float_end                         —          —         Y
    downsample_target_end                       params.downsample.target_end                        —          —         Y
    num_pyramid_layers                          params.num_pyramid_layers                           —          —         Y
    ==========================================  ================================================  =======  =========  =======

    Passing a kwarg not in the Y column for a given function raises TypeError.

    [1] pyramid_register defaults correspondences_flag_threshold to 0.999
        (MATLAB convention) when not explicitly passed. rigid/nonrigid use the
        C++ struct default.
    [2] pyramid_register auto-populates transform_num_viscous_iterations_start
        and transform_num_elastic_iterations_start to num_iterations when not
        explicitly passed. Explicit pass detected via kwargs.keys() — do not
        refactor to .get().

    Returns
    -------
    RigidRegResult
    """
    feat_float, feat_target, faces_float, faces_target, flags_float, flags_target = (
        _prepare_arrays(
            floating,
            target,
            floating_features=floating_features,
            target_features=target_features,
            floating_faces=floating_faces,
            target_faces=target_faces,
            floating_normals_override=normals,
            force_recompute_normals=compute_normals_flag,
            floating_flags=floating_flags,
            target_flags=target_flags,
        )
    )

    params = RigidParams()
    params = _apply_rigid_kwargs(params, kwargs)

    raw = _rigid_registration(
        feat_float,
        feat_target,
        faces_float,
        faces_target,
        flags_float,
        flags_target,
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
    rigid_params=None,
    **kwargs: Unpack[NonrigidKwargs],
) -> NonrigidRegResult:
    """Run nonrigid (viscoelastic) mesh registration.

    Accepts the same two call patterns and **kwargs as rigid_register; see
    that docstring for the full parameter table.

    Parameters
    ----------
    rigid_params:
        Optional dict of kwargs passed to rigid_register() for rigid pre-alignment.
        Pass ``{}`` to run rigid pre-alignment with defaults, or a dict with overrides
        such as ``{'num_iterations': 40}``.  ``None`` (default) skips pre-alignment.

    Returns
    -------
    NonrigidRegResult
    """
    feat_float, feat_target, faces_float, faces_target, flags_float, flags_target = (
        _prepare_arrays(
            floating,
            target,
            floating_features=floating_features,
            target_features=target_features,
            floating_faces=floating_faces,
            target_faces=target_faces,
            floating_normals_override=normals,
            force_recompute_normals=compute_normals_flag,
            floating_flags=floating_flags,
            target_flags=target_flags,
        )
    )

    if rigid_params is not None:
        if not isinstance(rigid_params, dict):
            raise TypeError(
                f"rigid_params must be a dict or None, got {type(rigid_params).__name__}"
            )
        _rigid_kw = rigid_params
        rigid_result = rigid_register(
            floating_features=feat_float,
            target_features=feat_target,
            floating_faces=faces_float,
            target_faces=faces_target,
            floating_flags=flags_float,
            target_flags=flags_target,
            **_rigid_kw,
        )
        feat_float = rigid_result.aligned_features

    params = NonrigidParams()
    params = _apply_nonrigid_kwargs(params, kwargs)

    raw = _nonrigid_registration(
        feat_float,
        feat_target,
        faces_float,
        faces_target,
        flags_float,
        flags_target,
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
    rigid_params=None,
    **kwargs: Unpack[PyramidKwargs],
) -> PyramidRegResult:
    """Run pyramid (multi-resolution) nonrigid mesh registration.

    Accepts the same two call patterns and **kwargs as rigid_register; see
    that docstring for the full parameter table.

    Special behaviour: viscous/elastic start iterations are auto-set to
    num_iterations when not explicitly passed (matching MATLAB convention).

    Parameters
    ----------
    rigid_params:
        Optional dict of kwargs passed to rigid_register() for rigid pre-alignment.
        Pass ``{}`` to run rigid pre-alignment with defaults, or a dict with overrides
        such as ``{'num_iterations': 40}``.  ``None`` (default) skips pre-alignment.

    Returns
    -------
    PyramidRegResult
    """
    feat_float, feat_target, faces_float, faces_target, flags_float, flags_target = (
        _prepare_arrays(
            floating,
            target,
            floating_features=floating_features,
            target_features=target_features,
            floating_faces=floating_faces,
            target_faces=target_faces,
            floating_normals_override=normals,
            force_recompute_normals=compute_normals_flag,
            floating_flags=floating_flags,
            target_flags=target_flags,
        )
    )

    if rigid_params is not None:
        if not isinstance(rigid_params, dict):
            raise TypeError(
                f"rigid_params must be a dict or None, got {type(rigid_params).__name__}"
            )
        _rigid_kw = rigid_params
        rigid_result = rigid_register(
            floating_features=feat_float,
            target_features=feat_target,
            floating_faces=faces_float,
            target_faces=faces_target,
            floating_flags=flags_float,
            target_flags=flags_target,
            **_rigid_kw,
        )
        feat_float = rigid_result.aligned_features

    # Pass the set of explicitly provided kwargs so _apply_pyramid_kwargs
    # knows which viscous/elastic start values to auto-populate
    params = PyramidParams()
    # Set pyramid-specific flag_threshold default (0.999 per MATLAB demo).
    # This is done here rather than in the C++ struct to keep PyramidParams
    # aggregate (no user-defined constructor), per ADR D6.
    if "correspondences_flag_threshold" not in kwargs:
        params.correspondences.flag_threshold = 0.999
    params = _apply_pyramid_kwargs(params, kwargs, explicit_kwargs=set(kwargs.keys()))

    raw = _pyramid_registration(
        feat_float,
        feat_target,
        faces_float,
        faces_target,
        flags_float,
        flags_target,
        params,
    )
    return PyramidRegResult(
        aligned_features=np.asarray(raw.aligned_features, dtype=np.float32),
        final_inlier_weights=np.asarray(raw.final_inlier_weights, dtype=np.float32),
        displacement_field=np.asarray(raw.displacement_field, dtype=np.float32),
        per_layer_iterations=list(raw.per_layer_iterations),
    )
