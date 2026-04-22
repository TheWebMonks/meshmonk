"""Type definitions for meshmonk register function **kwargs."""

from typing import TypedDict


class RigidKwargs(TypedDict, total=False):
    """TypedDict for **kwargs accepted by rigid_register().

    All fields are optional (total=False).  Fields match _RIGID_KNOWN_KWARGS
    in meshmonk/__init__.py.
    """

    correspondences_symmetric: bool
    correspondences_num_neighbours: int
    correspondences_flag_threshold: float
    correspondences_equalize_push_pull: bool
    inlier_kappa: float
    inlier_use_orientation: bool
    num_iterations: int
    use_scaling: bool


class NonrigidKwargs(TypedDict, total=False):
    """TypedDict for **kwargs accepted by nonrigid_register().

    All fields are optional (total=False).  Fields match _NONRIGID_KNOWN_KWARGS
    in meshmonk/__init__.py.
    """

    correspondences_symmetric: bool
    correspondences_num_neighbours: int
    correspondences_flag_threshold: float
    correspondences_equalize_push_pull: bool
    inlier_kappa: float
    inlier_use_orientation: bool
    transform_sigma: float
    transform_num_viscous_iterations_start: int
    transform_num_viscous_iterations_end: int
    transform_num_elastic_iterations_start: int
    transform_num_elastic_iterations_end: int
    num_iterations: int


class PyramidKwargs(TypedDict, total=False):
    """TypedDict for **kwargs accepted by pyramid_register().

    All fields are optional (total=False).  Fields match _PYRAMID_KNOWN_KWARGS
    in meshmonk/__init__.py.
    """

    correspondences_symmetric: bool
    correspondences_num_neighbours: int
    correspondences_flag_threshold: float
    correspondences_equalize_push_pull: bool
    inlier_kappa: float
    inlier_use_orientation: bool
    transform_sigma: float
    transform_num_viscous_iterations_start: int
    transform_num_viscous_iterations_end: int
    transform_num_elastic_iterations_start: int
    transform_num_elastic_iterations_end: int
    downsample_float_start: float
    downsample_target_start: float
    downsample_float_end: float
    downsample_target_end: float
    num_iterations: int
    num_pyramid_layers: int
