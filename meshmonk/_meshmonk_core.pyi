"""MeshMonk C++ extension module (nanobind)"""

import enum
from typing import Annotated

import numpy
from numpy.typing import NDArray

class MeshMonkError(RuntimeError):
    pass

class RegistrationError(enum.Enum):
    DegenerateInput = 0

    InsufficientInliers = 3

    DecompositionFailed = 1

    NonConvergence = 2

class LogLevel(enum.Enum):
    Silent = 0

    Error = 1

    Warning = 2

    Info = 3

    Debug = 4

class CorrespondenceParams:
    def __init__(self) -> None: ...
    @property
    def symmetric(self) -> bool: ...
    @symmetric.setter
    def symmetric(self, arg: bool, /) -> None: ...
    @property
    def num_neighbours(self) -> int: ...
    @num_neighbours.setter
    def num_neighbours(self, arg: int, /) -> None: ...
    @property
    def flag_threshold(self) -> float: ...
    @flag_threshold.setter
    def flag_threshold(self, arg: float, /) -> None: ...
    @property
    def equalize_push_pull(self) -> bool: ...
    @equalize_push_pull.setter
    def equalize_push_pull(self, arg: bool, /) -> None: ...

class InlierParams:
    def __init__(self) -> None: ...
    @property
    def kappa(self) -> float: ...
    @kappa.setter
    def kappa(self, arg: float, /) -> None: ...
    @property
    def use_orientation(self) -> bool: ...
    @use_orientation.setter
    def use_orientation(self, arg: bool, /) -> None: ...

class ViscoElasticParams:
    def __init__(self) -> None: ...
    @property
    def sigma(self) -> float: ...
    @sigma.setter
    def sigma(self, arg: float, /) -> None: ...
    @property
    def num_viscous_iterations_start(self) -> int: ...
    @num_viscous_iterations_start.setter
    def num_viscous_iterations_start(self, arg: int, /) -> None: ...
    @property
    def num_viscous_iterations_end(self) -> int: ...
    @num_viscous_iterations_end.setter
    def num_viscous_iterations_end(self, arg: int, /) -> None: ...
    @property
    def num_elastic_iterations_start(self) -> int: ...
    @num_elastic_iterations_start.setter
    def num_elastic_iterations_start(self, arg: int, /) -> None: ...
    @property
    def num_elastic_iterations_end(self) -> int: ...
    @num_elastic_iterations_end.setter
    def num_elastic_iterations_end(self, arg: int, /) -> None: ...

class DownsampleSchedule:
    def __init__(self) -> None: ...
    @property
    def float_start(self) -> float: ...
    @float_start.setter
    def float_start(self, arg: float, /) -> None: ...
    @property
    def target_start(self) -> float: ...
    @target_start.setter
    def target_start(self, arg: float, /) -> None: ...
    @property
    def float_end(self) -> float: ...
    @float_end.setter
    def float_end(self, arg: float, /) -> None: ...
    @property
    def target_end(self) -> float: ...
    @target_end.setter
    def target_end(self, arg: float, /) -> None: ...

class RigidParams:
    def __init__(self) -> None: ...
    @property
    def correspondences(self) -> CorrespondenceParams: ...
    @correspondences.setter
    def correspondences(self, arg: CorrespondenceParams, /) -> None: ...
    @property
    def inliers(self) -> InlierParams: ...
    @inliers.setter
    def inliers(self, arg: InlierParams, /) -> None: ...
    @property
    def num_iterations(self) -> int: ...
    @num_iterations.setter
    def num_iterations(self, arg: int, /) -> None: ...
    @property
    def use_scaling(self) -> bool: ...
    @use_scaling.setter
    def use_scaling(self, arg: bool, /) -> None: ...

class NonrigidParams:
    def __init__(self) -> None: ...
    @property
    def correspondences(self) -> CorrespondenceParams: ...
    @correspondences.setter
    def correspondences(self, arg: CorrespondenceParams, /) -> None: ...
    @property
    def inliers(self) -> InlierParams: ...
    @inliers.setter
    def inliers(self, arg: InlierParams, /) -> None: ...
    @property
    def transform(self) -> ViscoElasticParams: ...
    @transform.setter
    def transform(self, arg: ViscoElasticParams, /) -> None: ...
    @property
    def num_iterations(self) -> int: ...
    @num_iterations.setter
    def num_iterations(self, arg: int, /) -> None: ...

class PyramidParams:
    def __init__(self) -> None: ...
    @property
    def correspondences(self) -> CorrespondenceParams: ...
    @correspondences.setter
    def correspondences(self, arg: CorrespondenceParams, /) -> None: ...
    @property
    def inliers(self) -> InlierParams: ...
    @inliers.setter
    def inliers(self, arg: InlierParams, /) -> None: ...
    @property
    def transform(self) -> ViscoElasticParams: ...
    @transform.setter
    def transform(self, arg: ViscoElasticParams, /) -> None: ...
    @property
    def downsample(self) -> DownsampleSchedule: ...
    @downsample.setter
    def downsample(self, arg: DownsampleSchedule, /) -> None: ...
    @property
    def num_iterations(self) -> int: ...
    @num_iterations.setter
    def num_iterations(self, arg: int, /) -> None: ...
    @property
    def num_pyramid_layers(self) -> int: ...
    @num_pyramid_layers.setter
    def num_pyramid_layers(self, arg: int, /) -> None: ...

class RigidTransform:
    def __init__(self) -> None: ...
    @property
    def matrix(
        self,
    ) -> Annotated[NDArray[numpy.float32], dict(shape=(4, 4), order="F")]: ...
    def compose(self, rhs: RigidTransform) -> RigidTransform: ...
    def apply(
        self,
        features: Annotated[NDArray[numpy.float32], dict(shape=(None, 6), order="F")],
    ) -> Annotated[NDArray[numpy.float32], dict(shape=(None, 6), order="F")]: ...
    def inverse(self) -> RigidTransform: ...

class RigidResult:
    def __init__(self) -> None: ...
    @property
    def aligned_features(
        self,
    ) -> Annotated[NDArray[numpy.float32], dict(shape=(None, 6), order="F")]: ...
    @property
    def transform(self) -> RigidTransform: ...
    @property
    def iterations_run(self) -> int: ...

class NonrigidResult:
    def __init__(self) -> None: ...
    @property
    def aligned_features(
        self,
    ) -> Annotated[NDArray[numpy.float32], dict(shape=(None, 6), order="F")]: ...
    @property
    def final_inlier_weights(
        self,
    ) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order="C")]: ...
    @property
    def displacement_field(
        self,
    ) -> Annotated[NDArray[numpy.float32], dict(shape=(None, 3), order="F")]: ...
    @property
    def iterations_run(self) -> int: ...

class PyramidResult:
    def __init__(self) -> None: ...
    @property
    def aligned_features(
        self,
    ) -> Annotated[NDArray[numpy.float32], dict(shape=(None, 6), order="F")]: ...
    @property
    def final_inlier_weights(
        self,
    ) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order="C")]: ...
    @property
    def displacement_field(
        self,
    ) -> Annotated[NDArray[numpy.float32], dict(shape=(None, 3), order="F")]: ...
    @property
    def per_layer_iterations(self) -> list[int]: ...

def rigid_registration(
    floating_features: Annotated[
        NDArray[numpy.float32], dict(shape=(None, 6), writable=False)
    ],
    target_features: Annotated[
        NDArray[numpy.float32], dict(shape=(None, 6), writable=False)
    ],
    floating_faces: Annotated[
        NDArray[numpy.int32], dict(shape=(None, 3), writable=False)
    ],
    target_faces: Annotated[
        NDArray[numpy.int32], dict(shape=(None, 3), writable=False)
    ],
    floating_flags: Annotated[
        NDArray[numpy.float32], dict(shape=(None,), writable=False)
    ],
    target_flags: Annotated[
        NDArray[numpy.float32], dict(shape=(None,), writable=False)
    ],
    params: RigidParams = ...,
) -> RigidResult: ...
def nonrigid_registration(
    floating_features: Annotated[
        NDArray[numpy.float32], dict(shape=(None, 6), writable=False)
    ],
    target_features: Annotated[
        NDArray[numpy.float32], dict(shape=(None, 6), writable=False)
    ],
    floating_faces: Annotated[
        NDArray[numpy.int32], dict(shape=(None, 3), writable=False)
    ],
    target_faces: Annotated[
        NDArray[numpy.int32], dict(shape=(None, 3), writable=False)
    ],
    floating_flags: Annotated[
        NDArray[numpy.float32], dict(shape=(None,), writable=False)
    ],
    target_flags: Annotated[
        NDArray[numpy.float32], dict(shape=(None,), writable=False)
    ],
    params: NonrigidParams = ...,
) -> NonrigidResult: ...
def pyramid_registration(
    floating_features: Annotated[
        NDArray[numpy.float32], dict(shape=(None, 6), writable=False)
    ],
    target_features: Annotated[
        NDArray[numpy.float32], dict(shape=(None, 6), writable=False)
    ],
    floating_faces: Annotated[
        NDArray[numpy.int32], dict(shape=(None, 3), writable=False)
    ],
    target_faces: Annotated[
        NDArray[numpy.int32], dict(shape=(None, 3), writable=False)
    ],
    floating_flags: Annotated[
        NDArray[numpy.float32], dict(shape=(None,), writable=False)
    ],
    target_flags: Annotated[
        NDArray[numpy.float32], dict(shape=(None,), writable=False)
    ],
    params: PyramidParams = ...,
) -> PyramidResult: ...
def compute_correspondences(
    floating_features: Annotated[
        NDArray[numpy.float32], dict(shape=(None, 6), writable=False)
    ],
    target_features: Annotated[
        NDArray[numpy.float32], dict(shape=(None, 6), writable=False)
    ],
    floating_flags: Annotated[
        NDArray[numpy.float32], dict(shape=(None,), writable=False)
    ],
    target_flags: Annotated[
        NDArray[numpy.float32], dict(shape=(None,), writable=False)
    ],
    symmetric: bool = True,
    num_neighbours: int = 3,
    flag_threshold: float = 0.8999999761581421,
    equalize_push_pull: bool = False,
) -> tuple[
    Annotated[NDArray[numpy.float32], dict(shape=(None, 6), order="F")],
    Annotated[NDArray[numpy.float32], dict(shape=(None,), order="C")],
]: ...
def compute_inlier_weights(
    floating_features: Annotated[
        NDArray[numpy.float32], dict(shape=(None, 6), writable=False)
    ],
    corresponding_features: Annotated[
        NDArray[numpy.float32], dict(shape=(None, 6), writable=False)
    ],
    corresponding_flags: Annotated[
        NDArray[numpy.float32], dict(shape=(None,), writable=False)
    ],
    kappa: float = 12.0,
    use_orientation: bool = True,
) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order="C")]: ...
def compute_rigid_transform(
    floating_features: Annotated[
        NDArray[numpy.float32], dict(shape=(None, 6), writable=False)
    ],
    corresponding_features: Annotated[
        NDArray[numpy.float32], dict(shape=(None, 6), writable=False)
    ],
    weights: Annotated[NDArray[numpy.float32], dict(shape=(None,), writable=False)],
    use_scaling: bool = False,
) -> RigidTransform: ...
def compute_nonrigid_transform(
    floating_features: Annotated[
        NDArray[numpy.float32], dict(shape=(None, 6), writable=False)
    ],
    corresponding_features: Annotated[
        NDArray[numpy.float32], dict(shape=(None, 6), writable=False)
    ],
    floating_faces: Annotated[
        NDArray[numpy.int32], dict(shape=(None, 3), writable=False)
    ],
    floating_flags: Annotated[
        NDArray[numpy.float32], dict(shape=(None,), writable=False)
    ],
    weights: Annotated[NDArray[numpy.float32], dict(shape=(None,), writable=False)],
    num_smoothing_neighbours: int = 10,
    sigma: float = 3.0,
    num_viscous_iterations: int = 50,
    num_elastic_iterations: int = 50,
) -> Annotated[NDArray[numpy.float32], dict(shape=(None, 6), order="F")]: ...
def downsample_mesh(
    features: Annotated[NDArray[numpy.float32], dict(shape=(None, 6), writable=False)],
    faces: Annotated[NDArray[numpy.int32], dict(shape=(None, 3), writable=False)],
    flags: Annotated[NDArray[numpy.float32], dict(shape=(None,), writable=False)],
    downsample_ratio: float = 0.800000011920929,
) -> tuple[
    Annotated[NDArray[numpy.float32], dict(shape=(None, 6), order="F")],
    Annotated[NDArray[numpy.int32], dict(shape=(None, 3), order="F")],
    Annotated[NDArray[numpy.float32], dict(shape=(None,), order="C")],
    Annotated[NDArray[numpy.int32], dict(shape=(None,), order="C")],
]: ...
def scale_shift_mesh(
    previous_features: Annotated[
        NDArray[numpy.float32], dict(shape=(None, 6), writable=False)
    ],
    previous_indices: Annotated[
        NDArray[numpy.int32], dict(shape=(None,), writable=False)
    ],
    new_indices: Annotated[NDArray[numpy.int32], dict(shape=(None,), writable=False)],
) -> Annotated[NDArray[numpy.float32], dict(shape=(None, 6), order="F")]: ...
def compute_normals(
    positions: Annotated[NDArray[numpy.float32], dict(shape=(None, 3), writable=False)],
    faces: Annotated[NDArray[numpy.int32], dict(shape=(None, 3), writable=False)],
) -> Annotated[NDArray[numpy.float32], dict(shape=(None, 3), order="F")]: ...
def set_log_level(level: LogLevel) -> None: ...
