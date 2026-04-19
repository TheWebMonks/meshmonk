// Copyright 2026 jsnyde0
// SPDX-License-Identifier: Apache-2.0
//
// meshmonk.cpp — implements the public API declared in meshmonk/meshmonk.hpp.
// All internal classes live in namespace registration (load-bearing code, do NOT modify).
// This file bridges the new value-based public API to the old pointer-based internal API.

#include "meshmonk/meshmonk.hpp"

// Internal class headers (in library/src/ — included via target_include_directories)
#include "RigidRegistration.hpp"
#include "NonrigidRegistration.hpp"
#include "PyramidNonrigidRegistration.hpp"
#include "CorrespondenceFilter.hpp"
#include "SymmetricCorrespondenceFilter.hpp"
#include "BaseCorrespondenceFilter.hpp"
#include "InlierDetector.hpp"
#include "RigidTransformer.hpp"
#include "ViscoElasticTransformer.hpp"
#include "Downsampler.hpp"
#include "ScaleShifter.hpp"
#include "helper_functions.hpp"
#include "meshmonk/global.hpp"

#include <cmath>
#include <memory>

// Compile-time guard: both headers must agree on the feature count.
static_assert(meshmonk::NUM_FEATURES == registration::NUM_FEATURES,
              "NUM_FEATURES mismatch between meshmonk::types.hpp and global.hpp");

namespace meshmonk {

//=============================================================================
// Internal helpers
//=============================================================================

namespace {

// Run a single correspondence + inlier pass on already-aligned features.
// Used for final_inlier_weights in NonrigidResult and PyramidResult.
VecDynFloat run_final_inlier_pass(
    FeatureMat&        aligned_features,
    const FeatureMat&  target_copy,
    const VecDynFloat& floating_flags_copy,
    const VecDynFloat& target_flags_copy,
    const CorrespondenceParams& corr_params,
    const InlierParams&         inlier_params)
{
    const Eigen::Index N = aligned_features.rows();

    // Allocate output buffers
    FeatureMat  corresponding_features = FeatureMat::Zero(N, registration::NUM_FEATURES);
    VecDynFloat corresponding_flags    = VecDynFloat::Zero(N);
    VecDynFloat final_weights          = VecDynFloat::Ones(N);

    // Correspondence filter — branch on symmetric
    std::unique_ptr<registration::BaseCorrespondenceFilter> corr_filter;
    if (corr_params.symmetric) {
        auto f = std::make_unique<registration::SymmetricCorrespondenceFilter>();
        f->set_parameters(
            (size_t)corr_params.num_neighbours,
            corr_params.flag_threshold,
            corr_params.equalize_push_pull);
        corr_filter = std::move(f);
    } else {
        auto f = std::make_unique<registration::CorrespondenceFilter>();
        f->set_parameters(
            (size_t)corr_params.num_neighbours,
            corr_params.flag_threshold);
        corr_filter = std::move(f);
    }
    corr_filter->set_floating_input(&aligned_features, &floating_flags_copy);
    corr_filter->set_target_input(&target_copy, &target_flags_copy);
    corr_filter->set_output(&corresponding_features, &corresponding_flags);
    corr_filter->update();

    // Inlier detector — set_output() before update(); there is NO get_output()
    registration::InlierDetector inlier_det;
    inlier_det.set_input(&aligned_features, &corresponding_features, &corresponding_flags);
    inlier_det.set_output(&final_weights);   // output written in-place to final_weights
    inlier_det.set_parameters(inlier_params.kappa, inlier_params.use_orientation);
    inlier_det.update();

    return final_weights;
}

} // anonymous namespace


//=============================================================================
// High-level registration pipelines
//=============================================================================

[[nodiscard]] tl::expected<RigidResult, RegistrationError>
rigid_registration(
    Eigen::Ref<const FeatureMat>  floating_features,
    Eigen::Ref<const FeatureMat>  target_features,
    Eigen::Ref<const FacesMat>    floating_faces,
    Eigen::Ref<const FacesMat>    target_faces,
    Eigen::Ref<const VecDynFloat> floating_flags,
    Eigen::Ref<const VecDynFloat> target_flags,
    const RigidParams&            params)
{
    //-------------------------------------------------------------------------
    // DegenerateInput checks (before any algorithmic code)
    //-------------------------------------------------------------------------
    if (floating_features.rows() == 0 || target_features.rows() == 0)
        return tl::unexpected{RegistrationError::DegenerateInput};
    if (floating_features.rows() != floating_flags.size())
        return tl::unexpected{RegistrationError::DegenerateInput};
    if (target_features.rows() != target_flags.size())
        return tl::unexpected{RegistrationError::DegenerateInput};
    if (floating_faces.cols() != 3 || target_faces.cols() != 3)
        return tl::unexpected{RegistrationError::DegenerateInput};
    if (floating_faces.rows() == 0 || target_faces.rows() == 0)
        return tl::unexpected{RegistrationError::DegenerateInput};
    if (floating_flags.maxCoeff() == 0.0f)
        return tl::unexpected{RegistrationError::DegenerateInput};
    if (target_flags.maxCoeff() == 0.0f)
        return tl::unexpected{RegistrationError::DegenerateInput};
    // Zero-range bounding box check — degenerate only if ALL three axes have zero range
    if (floating_features.col(0).maxCoeff() == floating_features.col(0).minCoeff() &&
        floating_features.col(1).maxCoeff() == floating_features.col(1).minCoeff() &&
        floating_features.col(2).maxCoeff() == floating_features.col(2).minCoeff())
        return tl::unexpected{RegistrationError::DegenerateInput};
    // Target bounding box degeneracy check
    if (target_features.col(0).maxCoeff() == target_features.col(0).minCoeff() &&
        target_features.col(1).maxCoeff() == target_features.col(1).minCoeff() &&
        target_features.col(2).maxCoeff() == target_features.col(2).minCoeff())
        return tl::unexpected{RegistrationError::DegenerateInput};

    //-------------------------------------------------------------------------
    // Make mutable copies (set_input takes raw pointers to non-const)
    //-------------------------------------------------------------------------
    FeatureMat  floating_copy      = floating_features;
    FeatureMat  target_copy        = target_features;
    VecDynFloat floating_flags_copy = floating_flags;
    VecDynFloat target_flags_copy   = target_flags;

    //-------------------------------------------------------------------------
    // Run rigid registration — ALWAYS use set_parameters() explicitly
    //-------------------------------------------------------------------------
    registration::RigidRegistration registrator;
    registrator.set_input(&floating_copy, &target_copy, &floating_flags_copy, &target_flags_copy);
    registrator.set_parameters(
        params.correspondences.symmetric,
        (size_t)params.correspondences.num_neighbours,
        params.correspondences.flag_threshold,
        params.correspondences.equalize_push_pull,
        params.inliers.kappa,
        params.inliers.use_orientation,
        (size_t)params.num_iterations,
        params.use_scaling);
    registrator.update();

    //-------------------------------------------------------------------------
    // DecompositionFailed — detect post-hoc via allFinite()
    // (RigidTransformer logs to std::cerr on failure but doesn't throw)
    //-------------------------------------------------------------------------
    Eigen::Matrix4f mat = registrator.get_transformation();
    if (!mat.allFinite())
        return tl::unexpected{RegistrationError::DecompositionFailed};

    // TODO (v0.4): add convergence criterion (NonConvergence error)

    //-------------------------------------------------------------------------
    // InsufficientInliers check (consistent with nonrigid/pyramid)
    //-------------------------------------------------------------------------
    VecDynFloat final_weights = run_final_inlier_pass(
        floating_copy, target_copy, floating_flags_copy, target_flags_copy,
        params.correspondences, params.inliers);
    if ((final_weights.array() > 0.0f).count() < 4)
        return tl::unexpected{RegistrationError::InsufficientInliers};

    RigidResult result;
    result.aligned_features = floating_copy;
    result.transform        = RigidTransform{mat};
    result.iterations_run   = params.num_iterations;
    return result;
}


[[nodiscard]] tl::expected<NonrigidResult, RegistrationError>
nonrigid_registration(
    Eigen::Ref<const FeatureMat>  floating_features,
    Eigen::Ref<const FeatureMat>  target_features,
    Eigen::Ref<const FacesMat>    floating_faces,
    Eigen::Ref<const FacesMat>    target_faces,
    Eigen::Ref<const VecDynFloat> floating_flags,
    Eigen::Ref<const VecDynFloat> target_flags,
    const NonrigidParams&         params)
{
    //-------------------------------------------------------------------------
    // DegenerateInput checks
    //-------------------------------------------------------------------------
    if (floating_features.rows() == 0 || target_features.rows() == 0)
        return tl::unexpected{RegistrationError::DegenerateInput};
    if (floating_features.rows() != floating_flags.size())
        return tl::unexpected{RegistrationError::DegenerateInput};
    if (target_features.rows() != target_flags.size())
        return tl::unexpected{RegistrationError::DegenerateInput};
    if (floating_faces.cols() != 3 || target_faces.cols() != 3)
        return tl::unexpected{RegistrationError::DegenerateInput};
    if (floating_faces.rows() == 0 || target_faces.rows() == 0)
        return tl::unexpected{RegistrationError::DegenerateInput};
    if (floating_flags.maxCoeff() == 0.0f)
        return tl::unexpected{RegistrationError::DegenerateInput};
    if (target_flags.maxCoeff() == 0.0f)
        return tl::unexpected{RegistrationError::DegenerateInput};
    if (floating_features.col(0).maxCoeff() == floating_features.col(0).minCoeff() &&
        floating_features.col(1).maxCoeff() == floating_features.col(1).minCoeff() &&
        floating_features.col(2).maxCoeff() == floating_features.col(2).minCoeff())
        return tl::unexpected{RegistrationError::DegenerateInput};
    // Target bounding box degeneracy check
    if (target_features.col(0).maxCoeff() == target_features.col(0).minCoeff() &&
        target_features.col(1).maxCoeff() == target_features.col(1).minCoeff() &&
        target_features.col(2).maxCoeff() == target_features.col(2).minCoeff())
        return tl::unexpected{RegistrationError::DegenerateInput};
    if (params.num_iterations < 2)
        return tl::unexpected{RegistrationError::DegenerateInput};

    //-------------------------------------------------------------------------
    // Make mutable copies
    //-------------------------------------------------------------------------
    FeatureMat  floating_copy        = floating_features;
    FeatureMat  target_copy          = target_features;
    FacesMat    floating_faces_copy  = floating_faces;
    VecDynFloat floating_flags_copy  = floating_flags;
    VecDynFloat target_flags_copy    = target_flags;

    // Save original positions BEFORE update() for displacement_field computation
    Vec3Mat original_positions = floating_copy.leftCols(3);

    //-------------------------------------------------------------------------
    // Run nonrigid registration — ALWAYS use set_parameters() explicitly
    //-------------------------------------------------------------------------
    registration::NonrigidRegistration registrator;
    registrator.set_input(&floating_copy, &target_copy,
                          &floating_faces_copy,
                          &floating_flags_copy, &target_flags_copy);
    registrator.set_parameters(
        params.correspondences.symmetric,
        (size_t)params.correspondences.num_neighbours,
        params.correspondences.flag_threshold,
        params.correspondences.equalize_push_pull,
        params.inliers.kappa,
        params.inliers.use_orientation,
        (size_t)params.num_iterations,
        params.transform.sigma,
        (size_t)params.transform.num_viscous_iterations_start,
        (size_t)params.transform.num_viscous_iterations_end,
        (size_t)params.transform.num_elastic_iterations_start,
        (size_t)params.transform.num_elastic_iterations_end);
    registrator.update();  // modifies floating_copy in place

    // TODO (v0.4): add convergence criterion (NonConvergence error)

    //-------------------------------------------------------------------------
    // Compute final_inlier_weights via one additional correspondence+inlier pass
    //-------------------------------------------------------------------------
    VecDynFloat final_weights = run_final_inlier_pass(
        floating_copy, target_copy, floating_flags_copy, target_flags_copy,
        params.correspondences, params.inliers);

    //-------------------------------------------------------------------------
    // InsufficientInliers check
    //-------------------------------------------------------------------------
    if ((final_weights.array() > 0.0f).count() < 4)
        return tl::unexpected{RegistrationError::InsufficientInliers};

    NonrigidResult result;
    result.aligned_features     = floating_copy;
    result.final_inlier_weights = final_weights;
    result.displacement_field   = floating_copy.leftCols(3) - original_positions;
    result.iterations_run       = params.num_iterations;
    return result;
}


[[nodiscard]] tl::expected<PyramidResult, RegistrationError>
pyramid_registration(
    Eigen::Ref<const FeatureMat>  floating_features,
    Eigen::Ref<const FeatureMat>  target_features,
    Eigen::Ref<const FacesMat>    floating_faces,
    Eigen::Ref<const FacesMat>    target_faces,
    Eigen::Ref<const VecDynFloat> floating_flags,
    Eigen::Ref<const VecDynFloat> target_flags,
    const PyramidParams&          params)
{
    //-------------------------------------------------------------------------
    // DegenerateInput checks
    //-------------------------------------------------------------------------
    if (floating_features.rows() == 0 || target_features.rows() == 0)
        return tl::unexpected{RegistrationError::DegenerateInput};
    if (floating_features.rows() != floating_flags.size())
        return tl::unexpected{RegistrationError::DegenerateInput};
    if (target_features.rows() != target_flags.size())
        return tl::unexpected{RegistrationError::DegenerateInput};
    if (floating_faces.cols() != 3 || target_faces.cols() != 3)
        return tl::unexpected{RegistrationError::DegenerateInput};
    if (floating_faces.rows() == 0 || target_faces.rows() == 0)
        return tl::unexpected{RegistrationError::DegenerateInput};
    if (floating_flags.maxCoeff() == 0.0f)
        return tl::unexpected{RegistrationError::DegenerateInput};
    if (target_flags.maxCoeff() == 0.0f)
        return tl::unexpected{RegistrationError::DegenerateInput};
    if (floating_features.col(0).maxCoeff() == floating_features.col(0).minCoeff() &&
        floating_features.col(1).maxCoeff() == floating_features.col(1).minCoeff() &&
        floating_features.col(2).maxCoeff() == floating_features.col(2).minCoeff())
        return tl::unexpected{RegistrationError::DegenerateInput};
    // Target bounding box degeneracy check
    if (target_features.col(0).maxCoeff() == target_features.col(0).minCoeff() &&
        target_features.col(1).maxCoeff() == target_features.col(1).minCoeff() &&
        target_features.col(2).maxCoeff() == target_features.col(2).minCoeff())
        return tl::unexpected{RegistrationError::DegenerateInput};
    if (params.num_pyramid_layers < 1)
        return tl::unexpected{RegistrationError::DegenerateInput};
    {
        int iters_per_layer = (int)std::round(
            (float)params.num_iterations / (float)params.num_pyramid_layers);
        if (iters_per_layer < 2)
            return tl::unexpected{RegistrationError::DegenerateInput};
    }

    //-------------------------------------------------------------------------
    // Make mutable copies
    //-------------------------------------------------------------------------
    FeatureMat  floating_copy        = floating_features;
    FeatureMat  target_copy          = target_features;
    FacesMat    floating_faces_copy  = floating_faces;
    FacesMat    target_faces_copy    = target_faces;
    VecDynFloat floating_flags_copy  = floating_flags;
    VecDynFloat target_flags_copy    = target_flags;

    // Save original positions BEFORE update() for displacement_field computation
    Vec3Mat original_positions = floating_copy.leftCols(3);

    //-------------------------------------------------------------------------
    // Run pyramid nonrigid registration — ALWAYS use set_parameters() explicitly
    //-------------------------------------------------------------------------
    registration::PyramidNonrigidRegistration registrator;
    registrator.set_input(floating_copy, target_copy,
                          floating_faces_copy, target_faces_copy,
                          floating_flags_copy, target_flags_copy);
    registrator.set_parameters(
        (size_t)params.num_iterations,
        (size_t)params.num_pyramid_layers,
        params.downsample.float_start,
        params.downsample.target_start,
        params.downsample.float_end,
        params.downsample.target_end,
        params.correspondences.symmetric,
        (size_t)params.correspondences.num_neighbours,
        params.correspondences.flag_threshold,
        params.correspondences.equalize_push_pull,
        params.inliers.kappa,
        params.inliers.use_orientation,
        params.transform.sigma,
        (size_t)params.transform.num_viscous_iterations_start,
        (size_t)params.transform.num_viscous_iterations_end,
        (size_t)params.transform.num_elastic_iterations_start,
        (size_t)params.transform.num_elastic_iterations_end);
    registrator.update();  // modifies floating_copy in place

    // TODO (v0.4): add convergence criterion (NonConvergence error)

    //-------------------------------------------------------------------------
    // Compute per_layer_iterations
    // (PyramidNonrigidRegistration does not expose per-layer counts)
    // NOTE: iters_per_layer is an approximation — rounding means
    //       iters_per_layer * num_pyramid_layers may differ from num_iterations
    //       (e.g., 100 iterations / 3 layers → 33 per layer → 33*3 = 99).
    //-------------------------------------------------------------------------
    int iters_per_layer = (int)std::round(
        (float)params.num_iterations / (float)params.num_pyramid_layers);

    //-------------------------------------------------------------------------
    // Compute final_inlier_weights via one additional correspondence+inlier pass
    // (run at finest layer, i.e. on floating_copy after update())
    //-------------------------------------------------------------------------
    VecDynFloat final_weights = run_final_inlier_pass(
        floating_copy, target_copy, floating_flags_copy, target_flags_copy,
        params.correspondences, params.inliers);

    //-------------------------------------------------------------------------
    // InsufficientInliers check
    //-------------------------------------------------------------------------
    if ((final_weights.array() > 0.0f).count() < 4)
        return tl::unexpected{RegistrationError::InsufficientInliers};

    PyramidResult result;
    result.aligned_features     = floating_copy;
    result.final_inlier_weights = final_weights;
    result.displacement_field   = floating_copy.leftCols(3) - original_positions;
    result.per_layer_iterations = std::vector<int>(params.num_pyramid_layers, iters_per_layer);
    return result;
}


//=============================================================================
// Low-level primitives
//=============================================================================

std::pair<FeatureMat, VecDynFloat>
compute_correspondences(
    Eigen::Ref<const FeatureMat>  floating_features,
    Eigen::Ref<const FeatureMat>  target_features,
    Eigen::Ref<const VecDynFloat> floating_flags,
    Eigen::Ref<const VecDynFloat> target_flags,
    bool  symmetric,
    int   num_neighbours,
    float flag_threshold,
    bool  equalize_push_pull)
{
    // Make mutable copies for internal API (which takes const pointers)
    FeatureMat  floating_copy       = floating_features;
    FeatureMat  target_copy         = target_features;
    VecDynFloat floating_flags_copy = floating_flags;
    VecDynFloat target_flags_copy   = target_flags;

    const Eigen::Index N = floating_copy.rows();
    FeatureMat  corresponding_features = FeatureMat::Zero(N, registration::NUM_FEATURES);
    VecDynFloat corresponding_flags    = VecDynFloat::Zero(N);

    std::unique_ptr<registration::BaseCorrespondenceFilter> filter;
    if (symmetric) {
        auto f = std::make_unique<registration::SymmetricCorrespondenceFilter>();
        f->set_parameters((size_t)num_neighbours, flag_threshold, equalize_push_pull);
        filter = std::move(f);
    } else {
        auto f = std::make_unique<registration::CorrespondenceFilter>();
        f->set_parameters((size_t)num_neighbours, flag_threshold);
        filter = std::move(f);
    }
    filter->set_floating_input(&floating_copy, &floating_flags_copy);
    filter->set_target_input(&target_copy, &target_flags_copy);
    filter->set_output(&corresponding_features, &corresponding_flags);
    filter->update();

    return {corresponding_features, corresponding_flags};
}


VecDynFloat
compute_inlier_weights(
    Eigen::Ref<const FeatureMat>  floating_features,
    Eigen::Ref<const FeatureMat>  corresponding_features,
    Eigen::Ref<const VecDynFloat> corresponding_flags,
    float kappa,
    bool  use_orientation)
{
    FeatureMat  floating_copy      = floating_features;
    FeatureMat  corr_features_copy = corresponding_features;
    VecDynFloat corr_flags_copy    = corresponding_flags;
    VecDynFloat inlier_weights     = VecDynFloat::Ones(floating_copy.rows());

    registration::InlierDetector detector;
    detector.set_input(&floating_copy, &corr_features_copy, &corr_flags_copy);
    detector.set_output(&inlier_weights);
    detector.set_parameters(kappa, use_orientation);
    detector.update();

    return inlier_weights;
}


RigidTransform
compute_rigid_transform(
    Eigen::Ref<const FeatureMat>  floating_features,
    Eigen::Ref<const FeatureMat>  corresponding_features,
    Eigen::Ref<const VecDynFloat> weights,
    bool use_scaling)
{
    // All 4 steps required — RigidTransformer uses set_output target as
    // BOTH input for centroid computation AND write target; omitting set_output = null crash.
    FeatureMat  floating_copy      = floating_features;
    FeatureMat  corresponding_copy = corresponding_features;
    VecDynFloat weights_copy       = weights;

    registration::RigidTransformer t;
    t.set_input(&corresponding_copy, &weights_copy);  // (1) correspondences + weights
    t.set_output(&floating_copy);                     // (2) REQUIRED: used as both centroid input and write target
    t.set_parameters(use_scaling);                    // (3)
    t.update();                                       // (4) transforms floating_copy in-place

    return RigidTransform{t.get_transformation()};
}


FeatureMat
compute_nonrigid_transform(
    Eigen::Ref<const FeatureMat>  floating_features,
    Eigen::Ref<const FeatureMat>  corresponding_features,
    Eigen::Ref<const FacesMat>    floating_faces,
    Eigen::Ref<const VecDynFloat> floating_flags,
    Eigen::Ref<const VecDynFloat> weights,
    int   num_smoothing_neighbours,
    float sigma,
    int   num_viscous_iterations,
    int   num_elastic_iterations)
{
    FeatureMat  floating_copy      = floating_features;
    FeatureMat  corr_features_copy = corresponding_features;
    FacesMat    faces_copy         = floating_faces;
    VecDynFloat flags_copy         = floating_flags;
    VecDynFloat weights_copy       = weights;

    registration::ViscoElasticTransformer transformer;
    transformer.set_input(&corr_features_copy, &weights_copy, &flags_copy, &faces_copy);
    transformer.set_output(&floating_copy);
    transformer.set_parameters(
        (size_t)num_smoothing_neighbours,
        sigma,
        (size_t)num_viscous_iterations,
        (size_t)num_elastic_iterations);
    transformer.update();

    return floating_copy;
}


std::tuple<FeatureMat, FacesMat, VecDynFloat, Eigen::VectorXi>
downsample_mesh(
    Eigen::Ref<const FeatureMat>  features,
    Eigen::Ref<const FacesMat>    faces,
    Eigen::Ref<const VecDynFloat> flags,
    float downsample_ratio)
{
    FeatureMat  features_copy = features;
    FacesMat    faces_copy    = faces;
    VecDynFloat flags_copy    = flags;

    FeatureMat       out_features;
    FacesMat         out_faces;
    VecDynFloat      out_flags;
    Eigen::VectorXi  out_indices;

    registration::Downsampler downsampler;
    downsampler.set_input(&features_copy, &faces_copy, &flags_copy);
    downsampler.set_output(out_features, out_faces, out_flags, out_indices);
    downsampler.set_parameters(downsample_ratio);
    downsampler.update();

    return {out_features, out_faces, out_flags, out_indices};
}


FeatureMat
scale_shift_mesh(
    Eigen::Ref<const FeatureMat>      previous_features,
    Eigen::Ref<const Eigen::VectorXi> previous_indices,
    Eigen::Ref<const Eigen::VectorXi> new_indices)
{
    FeatureMat       prev_features_copy = previous_features;
    Eigen::VectorXi  prev_indices_copy  = previous_indices;
    Eigen::VectorXi  new_indices_copy   = new_indices;

    FeatureMat out_features(new_indices_copy.size(), registration::NUM_FEATURES);

    registration::ScaleShifter scaler;
    scaler.set_input(prev_features_copy, prev_indices_copy, new_indices_copy);
    scaler.set_output(out_features);
    scaler.update();

    return out_features;
}


Vec3Mat
compute_normals(
    Eigen::Ref<const Vec3Mat>   positions,
    Eigen::Ref<const FacesMat>  faces)
{
    Vec3Mat  positions_copy = positions;
    FacesMat faces_copy     = faces;
    Vec3Mat  normals(positions_copy.rows(), 3);

    registration::update_normals_for_altered_positions(positions_copy, faces_copy, normals);

    return normals;
}

} // namespace meshmonk
