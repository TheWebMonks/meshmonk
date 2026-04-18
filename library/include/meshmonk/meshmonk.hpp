#ifndef MESHMONK_HPP
#define MESHMONK_HPP

#include <utility>
#include <tuple>
#include <tl/expected.hpp>
#include <stdexcept>
#include <string>

#include "types.hpp"
#include "transform.hpp"
#include "params.hpp"
#include "result.hpp"
#include "logger.hpp"

namespace meshmonk {

class MeshMonkError : public std::runtime_error {
public:
    explicit MeshMonkError(RegistrationError code, const std::string& msg)
        : std::runtime_error(msg), code_(code) {}
    RegistrationError code() const noexcept { return code_; }
private:
    RegistrationError code_;
};

//=============================================================================
// High-level registration pipelines
// All return tl::expected<Result, RegistrationError> — callers MUST check.
//=============================================================================

[[nodiscard]] tl::expected<RigidResult, RegistrationError>
rigid_registration(
    Eigen::Ref<const FeatureMat>  floating_features,
    Eigen::Ref<const FeatureMat>  target_features,
    Eigen::Ref<const FacesMat>    floating_faces,
    Eigen::Ref<const FacesMat>    target_faces,
    Eigen::Ref<const VecDynFloat> floating_flags,
    Eigen::Ref<const VecDynFloat> target_flags,
    const RigidParams&            params = {});

[[nodiscard]] tl::expected<NonrigidResult, RegistrationError>
nonrigid_registration(
    Eigen::Ref<const FeatureMat>  floating_features,
    Eigen::Ref<const FeatureMat>  target_features,
    Eigen::Ref<const FacesMat>    floating_faces,
    Eigen::Ref<const FacesMat>    target_faces,
    Eigen::Ref<const VecDynFloat> floating_flags,
    Eigen::Ref<const VecDynFloat> target_flags,
    const NonrigidParams&         params = {});

[[nodiscard]] tl::expected<PyramidResult, RegistrationError>
pyramid_registration(
    Eigen::Ref<const FeatureMat>  floating_features,
    Eigen::Ref<const FeatureMat>  target_features,
    Eigen::Ref<const FacesMat>    floating_faces,
    Eigen::Ref<const FacesMat>    target_faces,
    Eigen::Ref<const VecDynFloat> floating_flags,
    Eigen::Ref<const VecDynFloat> target_flags,
    const PyramidParams&          params = {});


//=============================================================================
// Low-level primitives — value-returning (no in-place mutation of inputs)
//=============================================================================

// Compute weighted correspondences in 6D feature space (positions + normals)
// Returns: {corresponding_features (N,6), corresponding_flags (N,)}
std::pair<FeatureMat, VecDynFloat>
compute_correspondences(
    Eigen::Ref<const FeatureMat>  floating_features,
    Eigen::Ref<const FeatureMat>  target_features,
    Eigen::Ref<const VecDynFloat> floating_flags,
    Eigen::Ref<const VecDynFloat> target_flags,
    bool  symmetric          = true,
    int   num_neighbours     = 3,
    float flag_threshold     = 0.9f,
    bool  equalize_push_pull = false);

// Compute inlier weights given features and their correspondences
// Returns: inlier_weights (N,)
VecDynFloat
compute_inlier_weights(
    Eigen::Ref<const FeatureMat>  floating_features,
    Eigen::Ref<const FeatureMat>  corresponding_features,
    Eigen::Ref<const VecDynFloat> corresponding_flags,
    float kappa           = 12.0f,
    bool  use_orientation = true);

// Compute rigid transform from features + correspondences + weights.
// IMPORTANT: returns RigidTransform (does NOT mutate floating_features in place).
// Implementation requires set_output() on RigidTransformer (REQUIRED or null pointer crash).
RigidTransform
compute_rigid_transform(
    Eigen::Ref<const FeatureMat>  floating_features,
    Eigen::Ref<const FeatureMat>  corresponding_features,
    Eigen::Ref<const VecDynFloat> weights,
    bool use_scaling = false);

// Apply nonrigid deformation.
// Returns: deformed_features (N,6)
FeatureMat
compute_nonrigid_transform(
    Eigen::Ref<const FeatureMat>  floating_features,
    Eigen::Ref<const FeatureMat>  corresponding_features,
    Eigen::Ref<const FacesMat>    floating_faces,
    Eigen::Ref<const VecDynFloat> floating_flags,
    Eigen::Ref<const VecDynFloat> weights,
    int   num_smoothing_neighbours = 10,
    float sigma                    = 3.0f,
    int   num_viscous_iterations   = 50,
    int   num_elastic_iterations   = 50);

// Downsample a mesh.
// Returns: {features, faces, flags, original_indices}
std::tuple<FeatureMat, FacesMat, VecDynFloat, Eigen::VectorXi>
downsample_mesh(
    Eigen::Ref<const FeatureMat>  features,
    Eigen::Ref<const FacesMat>    faces,
    Eigen::Ref<const VecDynFloat> flags,
    float downsample_ratio = 0.8f);

// Transfer floating mesh features between pyramid scales.
// Returns: interpolated features for the new (higher-density) scale
FeatureMat
scale_shift_mesh(
    Eigen::Ref<const FeatureMat>      previous_features,
    Eigen::Ref<const Eigen::VectorXi> previous_indices,
    Eigen::Ref<const Eigen::VectorXi> new_indices);

// Compute per-vertex normals from positions + face topology.
// Returns: normals (N,3)
Vec3Mat
compute_normals(
    Eigen::Ref<const Vec3Mat>   positions,
    Eigen::Ref<const FacesMat>  faces);

} // namespace meshmonk

#endif // MESHMONK_HPP
