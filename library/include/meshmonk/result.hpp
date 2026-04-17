#ifndef MESHMONK_RESULT_HPP
#define MESHMONK_RESULT_HPP

#include <vector>
#include "types.hpp"
#include "transform.hpp"

namespace meshmonk {

enum class RegistrationError {
    DegenerateInput,     // empty inputs, dim/shape mismatch, zero-range bbox, all-zero flags
    DecompositionFailed, // RigidTransformer eigendecomposition fails (Horn method / SelfAdjointEigenSolver)
    NonConvergence,      // v0.2+ only — NOT raised in v0.1
    InsufficientInliers, // fewer than 4 non-zero inlier weights after InlierDetector
};

struct RigidResult {
    FeatureMat     aligned_features;     // (N,6) positions + normals
    RigidTransform transform;            // 4x4 SE(3)
    int            iterations_run = 0;
};

struct NonrigidResult {
    FeatureMat  aligned_features;      // (N,6)
    VecDynFloat final_inlier_weights;  // (N,) — computed via a final correspondence+inlier pass after registration
    Vec3Mat     displacement_field;    // (N,3) per-vertex displacement from original floating positions
    int         iterations_run = 0;
};

struct PyramidResult {
    FeatureMat         aligned_features;     // (N,6) at finest layer
    VecDynFloat        final_inlier_weights; // (N,) — computed via a final pass at finest layer
    Vec3Mat            displacement_field;   // (N,3) at finest layer
    std::vector<int>   per_layer_iterations;
};

} // namespace meshmonk

#endif // MESHMONK_RESULT_HPP
