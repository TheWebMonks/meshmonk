#ifndef MESHMONK_PARAMS_HPP
#define MESHMONK_PARAMS_HPP

// All param structs use MATLAB demo values as defaults (NOT legacy C++ defaults).
// Field order is STABLE across minor versions (v0.x -> v0.x+1); reordering requires an ADR update.
// C++ callers use C++20 designated initializers; nanobind shim uses named nb::arg() bindings.

namespace meshmonk {

struct CorrespondenceParams {
    bool  symmetric          = true;
    int   num_neighbours     = 3;      // MATLAB: 3, NOT legacy C++ default of 5
    float flag_threshold     = 0.9f;   // MATLAB: 0.9, NOT legacy C++ default of 0.99
    bool  equalize_push_pull = false;
};

struct InlierParams {
    float kappa           = 12.0f;  // MATLAB: 12.0, NOT legacy C++ default of 4.0
    bool  use_orientation = true;
};

struct ViscoElasticParams {
    float sigma                       = 3.0f;
    int   num_viscous_iterations_start = 200;  // fallback for nonrigid; pyramid auto-sets to num_iterations at Python level
    int   num_viscous_iterations_end   = 1;
    int   num_elastic_iterations_start = 200;
    int   num_elastic_iterations_end   = 1;
};

struct DownsampleSchedule {
    float float_start  = 50.f;
    float target_start = 70.f;
    float float_end    = 0.f;
    float target_end   = 0.f;
};

struct RigidParams {
    CorrespondenceParams correspondences;
    InlierParams         inliers;
    int                  num_iterations = 80;  // MATLAB: 80, NOT legacy C++ default of 20
    bool                 use_scaling    = false;
};

struct NonrigidParams {
    CorrespondenceParams correspondences;
    InlierParams         inliers;
    ViscoElasticParams   transform;
    int                  num_iterations = 200;  // MATLAB: 200, NOT legacy C++ default of 60
};

struct PyramidParams {
    CorrespondenceParams correspondences;
    InlierParams         inliers;
    ViscoElasticParams   transform;
    DownsampleSchedule   downsample;
    int                  num_iterations    = 90;  // MATLAB: 90, NOT legacy C++ default of 60
    int                  num_pyramid_layers = 3;

    // Override flag_threshold to MATLAB pyramid default (0.999).
    // CorrespondenceParams defaults to 0.9 which is correct for rigid/nonrigid,
    // but pyramid uses 0.999 (matching MATLAB demo).
    PyramidParams() {
        correspondences.flag_threshold = 0.999f;
    }
};

} // namespace meshmonk

#endif // MESHMONK_PARAMS_HPP
