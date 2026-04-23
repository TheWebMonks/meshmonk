// Copyright 2026 jsnyde0
// SPDX-License-Identifier: Apache-2.0
//
// bindings.cpp — nanobind extension module _meshmonk_core
//
// Boundary-copy note:
//   numpy arrays are row-major by default; Eigen is column-major.
//   nanobind performs one implicit copy per call at the Python->C++ boundary
//   (~1.3 MB for a 54K-vertex mesh). This is expected and acceptable.
//   Do not add workarounds — profiling in v0.2 will determine if it matters.
//
// GIL release note:
//   nb::call_guard<nb::gil_scoped_release>() releases the GIL for the
//   FUNCTION BODY only. Argument parsing (including Eigen::Ref type casters)
//   happens BEFORE the GIL is released, so it is safe to accept Eigen::Ref
//   parameters directly. Do NOT use nb::cast inside the function body as
//   that would call into Python with GIL released (segfault).

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/vector.h>

#include <meshmonk/logger.hpp>
#include <meshmonk/meshmonk.hpp>
#include <meshmonk/params.hpp>
#include <meshmonk/profiling.hpp>
#include <meshmonk/result.hpp>
#include <meshmonk/transform.hpp>

#include <stdexcept>
#include <string>

namespace nb = nanobind;
using namespace nb::literals;

// ---------------------------------------------------------------------------
// Error-unwrap helper: convert tl::expected failure to MeshMonkError
// Called inside the function body (after GIL is released) — only throws,
// does not interact with Python runtime.
// ---------------------------------------------------------------------------
template <typename T>
T unwrap_expected(tl::expected<T, meshmonk::RegistrationError> result) {
  if (result.has_value()) {
    return std::move(result.value());
  }
  meshmonk::RegistrationError code = result.error();
  std::string msg;
  switch (code) {
  case meshmonk::RegistrationError::DegenerateInput:
    msg = "DegenerateInput: empty inputs, shape mismatch, zero-range bounding "
          "box, or all-zero inlier flags";
    break;
  case meshmonk::RegistrationError::InsufficientInliers:
    msg = "InsufficientInliers: fewer than 4 non-zero inlier weights after "
          "InlierDetector";
    break;
  case meshmonk::RegistrationError::DecompositionFailed:
    msg = "DecompositionFailed: RigidTransformer eigendecomposition (Horn "
          "method) failed to converge";
    break;
  case meshmonk::RegistrationError::NonConvergence:
    msg = "NonConvergence: registration did not converge within the allowed "
          "iterations";
    break;
  default:
    msg = "Unknown RegistrationError";
    break;
  }
  throw meshmonk::MeshMonkError(code, msg);
}

NB_MODULE(_meshmonk_core, m) {
  m.doc() = "MeshMonk C++ extension module (nanobind)";

  // -----------------------------------------------------------------------
  // Exception translation: MeshMonkError -> Python MeshMonkError (RuntimeError)
  // We create the Python exception type manually so we can attach _code as
  // an integer attribute, enabling programmatic dispatch on error type.
  // -----------------------------------------------------------------------
  // exception_new returns a new reference; steal it into the module attr
  // so the refcount stays at 1 (owned by the module dict).
  // s_MeshMonkError is a borrowed raw pointer valid for the module lifetime.
  static PyObject *s_MeshMonkError =
      nb::detail::exception_new(m.ptr(), "MeshMonkError", PyExc_RuntimeError);

  // Export MeshMonkError on the module (steals the new reference)
  m.attr("MeshMonkError") = nb::steal(nb::handle(s_MeshMonkError));

  // Register custom translator that sets _code attribute on exception
  nb::register_exception_translator(
      [](const std::exception_ptr &p, void *payload) {
        try {
          std::rethrow_exception(p);
        } catch (const meshmonk::MeshMonkError &e) {
          auto *type = static_cast<PyObject *>(payload);
          PyObject *args = Py_BuildValue("(s)", e.what());
          if (!args)
            return;
          PyObject *exc = PyObject_Call(type, args, nullptr);
          Py_DECREF(args);
          if (exc) {
            PyObject *code_int = PyLong_FromLong(static_cast<long>(e.code()));
            PyObject_SetAttrString(exc, "_code", code_int);
            Py_XDECREF(code_int);
            PyErr_SetObject(type, exc);
            Py_DECREF(exc);
          }
        }
      },
      static_cast<void *>(s_MeshMonkError));

  // -----------------------------------------------------------------------
  // RegistrationError enum
  // -----------------------------------------------------------------------
  nb::enum_<meshmonk::RegistrationError>(m, "RegistrationError")
      .value("DegenerateInput", meshmonk::RegistrationError::DegenerateInput)
      .value("InsufficientInliers",
             meshmonk::RegistrationError::InsufficientInliers)
      .value("DecompositionFailed",
             meshmonk::RegistrationError::DecompositionFailed)
      .value("NonConvergence", meshmonk::RegistrationError::NonConvergence);

  // -----------------------------------------------------------------------
  // LogLevel enum
  // -----------------------------------------------------------------------
  nb::enum_<meshmonk::LogLevel>(m, "LogLevel")
      .value("Silent", meshmonk::LogLevel::Silent)
      .value("Error", meshmonk::LogLevel::Error)
      .value("Warning", meshmonk::LogLevel::Warning)
      .value("Info", meshmonk::LogLevel::Info)
      .value("Debug", meshmonk::LogLevel::Debug);

  // -----------------------------------------------------------------------
  // Param structs
  // -----------------------------------------------------------------------
  nb::class_<meshmonk::CorrespondenceParams>(m, "CorrespondenceParams")
      .def(nb::init<>())
      .def_rw("symmetric", &meshmonk::CorrespondenceParams::symmetric)
      .def_rw("num_neighbours", &meshmonk::CorrespondenceParams::num_neighbours)
      .def_rw("flag_threshold", &meshmonk::CorrespondenceParams::flag_threshold)
      .def_rw("equalize_push_pull",
              &meshmonk::CorrespondenceParams::equalize_push_pull);

  nb::class_<meshmonk::InlierParams>(m, "InlierParams")
      .def(nb::init<>())
      .def_rw("kappa", &meshmonk::InlierParams::kappa)
      .def_rw("use_orientation", &meshmonk::InlierParams::use_orientation);

  nb::class_<meshmonk::ViscoElasticParams>(m, "ViscoElasticParams")
      .def(nb::init<>())
      .def_rw("sigma", &meshmonk::ViscoElasticParams::sigma)
      .def_rw("num_viscous_iterations_start",
              &meshmonk::ViscoElasticParams::num_viscous_iterations_start)
      .def_rw("num_viscous_iterations_end",
              &meshmonk::ViscoElasticParams::num_viscous_iterations_end)
      .def_rw("num_elastic_iterations_start",
              &meshmonk::ViscoElasticParams::num_elastic_iterations_start)
      .def_rw("num_elastic_iterations_end",
              &meshmonk::ViscoElasticParams::num_elastic_iterations_end);

  nb::class_<meshmonk::DownsampleSchedule>(m, "DownsampleSchedule")
      .def(nb::init<>())
      .def_rw("float_start", &meshmonk::DownsampleSchedule::float_start)
      .def_rw("target_start", &meshmonk::DownsampleSchedule::target_start)
      .def_rw("float_end", &meshmonk::DownsampleSchedule::float_end)
      .def_rw("target_end", &meshmonk::DownsampleSchedule::target_end);

  nb::class_<meshmonk::RigidParams>(m, "RigidParams")
      .def(nb::init<>())
      .def_rw("correspondences", &meshmonk::RigidParams::correspondences)
      .def_rw("inliers", &meshmonk::RigidParams::inliers)
      .def_rw("num_iterations", &meshmonk::RigidParams::num_iterations)
      .def_rw("use_scaling", &meshmonk::RigidParams::use_scaling);

  nb::class_<meshmonk::NonrigidParams>(m, "NonrigidParams")
      .def(nb::init<>())
      .def_rw("correspondences", &meshmonk::NonrigidParams::correspondences)
      .def_rw("inliers", &meshmonk::NonrigidParams::inliers)
      .def_rw("transform", &meshmonk::NonrigidParams::transform)
      .def_rw("num_iterations", &meshmonk::NonrigidParams::num_iterations);

  nb::class_<meshmonk::PyramidParams>(m, "PyramidParams")
      .def(nb::init<>())
      .def_rw("correspondences", &meshmonk::PyramidParams::correspondences)
      .def_rw("inliers", &meshmonk::PyramidParams::inliers)
      .def_rw("transform", &meshmonk::PyramidParams::transform)
      .def_rw("downsample", &meshmonk::PyramidParams::downsample)
      .def_rw("num_iterations", &meshmonk::PyramidParams::num_iterations)
      .def_rw("num_pyramid_layers",
              &meshmonk::PyramidParams::num_pyramid_layers);

  // -----------------------------------------------------------------------
  // RigidTransform
  // -----------------------------------------------------------------------
  nb::class_<meshmonk::RigidTransform>(m, "RigidTransform")
      .def(nb::init<>())
      .def_prop_ro("matrix",
                   [](const meshmonk::RigidTransform &t) { return t.matrix; })
      .def("compose", &meshmonk::RigidTransform::compose, "rhs"_a)
      .def("apply", &meshmonk::RigidTransform::apply, "features"_a)
      .def("inverse", &meshmonk::RigidTransform::inverse);

  // -----------------------------------------------------------------------
  // Result structs (read-only from Python)
  // -----------------------------------------------------------------------
  nb::class_<meshmonk::RigidResult>(m, "RigidResult")
      .def(nb::init<>())
      .def_ro("aligned_features", &meshmonk::RigidResult::aligned_features)
      .def_ro("transform", &meshmonk::RigidResult::transform)
      .def_ro("iterations_run", &meshmonk::RigidResult::iterations_run);

  nb::class_<meshmonk::NonrigidResult>(m, "NonrigidResult")
      .def(nb::init<>())
      .def_ro("aligned_features", &meshmonk::NonrigidResult::aligned_features)
      .def_ro("final_inlier_weights",
              &meshmonk::NonrigidResult::final_inlier_weights)
      .def_ro("displacement_field",
              &meshmonk::NonrigidResult::displacement_field,
              "Per-vertex displacement from original floating mesh position. "
              "Shape (N, 3), float32.")
      .def_ro("iterations_run", &meshmonk::NonrigidResult::iterations_run);

  nb::class_<meshmonk::PyramidResult>(m, "PyramidResult")
      .def(nb::init<>())
      .def_ro("aligned_features", &meshmonk::PyramidResult::aligned_features)
      .def_ro("final_inlier_weights",
              &meshmonk::PyramidResult::final_inlier_weights)
      .def_ro("displacement_field",
              &meshmonk::PyramidResult::displacement_field,
              "Per-vertex displacement from original floating mesh position. "
              "Shape (N, 3), float32.")
      .def_ro("per_layer_iterations",
              &meshmonk::PyramidResult::per_layer_iterations);

  // -----------------------------------------------------------------------
  // High-level registration pipelines
  //
  // IMPORTANT: nb::call_guard<nb::gil_scoped_release>() releases the GIL
  // for the FUNCTION BODY only. Argument parsing (Eigen::Ref type casters)
  // happens BEFORE the GIL is released, so direct Eigen::Ref parameters
  // are safe. The function body then runs without the GIL.
  //
  // GIL release prevents blocking all Python threads during 30-180s calls.
  // nanobind GC-pins numpy arrays automatically so data lifetime is safe.
  // -----------------------------------------------------------------------
  m.def(
      "rigid_registration",
      [](Eigen::Ref<const meshmonk::FeatureMat> floating_features,
         Eigen::Ref<const meshmonk::FeatureMat> target_features,
         Eigen::Ref<const meshmonk::FacesMat> floating_faces,
         Eigen::Ref<const meshmonk::FacesMat> target_faces,
         Eigen::Ref<const meshmonk::VecDynFloat> floating_flags,
         Eigen::Ref<const meshmonk::VecDynFloat> target_flags,
         const meshmonk::RigidParams &params) -> meshmonk::RigidResult {
        return unwrap_expected(meshmonk::rigid_registration(
            floating_features, target_features, floating_faces, target_faces,
            floating_flags, target_flags, params));
      },
      "floating_features"_a, "target_features"_a, "floating_faces"_a,
      "target_faces"_a, "floating_flags"_a, "target_flags"_a,
      "params"_a = meshmonk::RigidParams{},
      nb::call_guard<nb::gil_scoped_release>());

  m.def(
      "nonrigid_registration",
      [](Eigen::Ref<const meshmonk::FeatureMat> floating_features,
         Eigen::Ref<const meshmonk::FeatureMat> target_features,
         Eigen::Ref<const meshmonk::FacesMat> floating_faces,
         Eigen::Ref<const meshmonk::FacesMat> target_faces,
         Eigen::Ref<const meshmonk::VecDynFloat> floating_flags,
         Eigen::Ref<const meshmonk::VecDynFloat> target_flags,
         const meshmonk::NonrigidParams &params) -> meshmonk::NonrigidResult {
        return unwrap_expected(meshmonk::nonrigid_registration(
            floating_features, target_features, floating_faces, target_faces,
            floating_flags, target_flags, params));
      },
      "floating_features"_a, "target_features"_a, "floating_faces"_a,
      "target_faces"_a, "floating_flags"_a, "target_flags"_a,
      "params"_a = meshmonk::NonrigidParams{},
      nb::call_guard<nb::gil_scoped_release>());

  m.def(
      "pyramid_registration",
      [](Eigen::Ref<const meshmonk::FeatureMat> floating_features,
         Eigen::Ref<const meshmonk::FeatureMat> target_features,
         Eigen::Ref<const meshmonk::FacesMat> floating_faces,
         Eigen::Ref<const meshmonk::FacesMat> target_faces,
         Eigen::Ref<const meshmonk::VecDynFloat> floating_flags,
         Eigen::Ref<const meshmonk::VecDynFloat> target_flags,
         const meshmonk::PyramidParams &params) -> meshmonk::PyramidResult {
        return unwrap_expected(meshmonk::pyramid_registration(
            floating_features, target_features, floating_faces, target_faces,
            floating_flags, target_flags, params));
      },
      "floating_features"_a, "target_features"_a, "floating_faces"_a,
      "target_faces"_a, "floating_flags"_a, "target_flags"_a,
      "params"_a = meshmonk::PyramidParams{},
      nb::call_guard<nb::gil_scoped_release>());

  // -----------------------------------------------------------------------
  // Low-level primitives (GIL released — safe for the same reason as above)
  // -----------------------------------------------------------------------
  m.def(
      "compute_correspondences",
      [](Eigen::Ref<const meshmonk::FeatureMat> floating_features,
         Eigen::Ref<const meshmonk::FeatureMat> target_features,
         Eigen::Ref<const meshmonk::VecDynFloat> floating_flags,
         Eigen::Ref<const meshmonk::VecDynFloat> target_flags, bool symmetric,
         int num_neighbours, float flag_threshold, bool equalize_push_pull)
          -> std::pair<meshmonk::FeatureMat, meshmonk::VecDynFloat> {
        return meshmonk::compute_correspondences(
            floating_features, target_features, floating_flags, target_flags,
            symmetric, num_neighbours, flag_threshold, equalize_push_pull);
      },
      "floating_features"_a, "target_features"_a, "floating_flags"_a,
      "target_flags"_a, "symmetric"_a = true, "num_neighbours"_a = 3,
      "flag_threshold"_a = 0.9f, "equalize_push_pull"_a = false,
      nb::call_guard<nb::gil_scoped_release>());

  m.def(
      "compute_inlier_weights",
      [](Eigen::Ref<const meshmonk::FeatureMat> floating_features,
         Eigen::Ref<const meshmonk::FeatureMat> corresponding_features,
         Eigen::Ref<const meshmonk::VecDynFloat> corresponding_flags,
         float kappa, bool use_orientation) -> meshmonk::VecDynFloat {
        return meshmonk::compute_inlier_weights(
            floating_features, corresponding_features, corresponding_flags,
            kappa, use_orientation);
      },
      "floating_features"_a, "corresponding_features"_a,
      "corresponding_flags"_a, "kappa"_a = 12.0f, "use_orientation"_a = true,
      nb::call_guard<nb::gil_scoped_release>());

  m.def(
      "compute_rigid_transform",
      [](Eigen::Ref<const meshmonk::FeatureMat> floating_features,
         Eigen::Ref<const meshmonk::FeatureMat> corresponding_features,
         Eigen::Ref<const meshmonk::VecDynFloat> weights,
         bool use_scaling) -> meshmonk::RigidTransform {
        return meshmonk::compute_rigid_transform(
            floating_features, corresponding_features, weights, use_scaling);
      },
      "floating_features"_a, "corresponding_features"_a, "weights"_a,
      "use_scaling"_a = false, nb::call_guard<nb::gil_scoped_release>());

  m.def(
      "compute_nonrigid_transform",
      [](Eigen::Ref<const meshmonk::FeatureMat> floating_features,
         Eigen::Ref<const meshmonk::FeatureMat> corresponding_features,
         Eigen::Ref<const meshmonk::FacesMat> floating_faces,
         Eigen::Ref<const meshmonk::VecDynFloat> floating_flags,
         Eigen::Ref<const meshmonk::VecDynFloat> weights,
         int num_smoothing_neighbours, float sigma, int num_viscous_iterations,
         int num_elastic_iterations) -> meshmonk::FeatureMat {
        return meshmonk::compute_nonrigid_transform(
            floating_features, corresponding_features, floating_faces,
            floating_flags, weights, num_smoothing_neighbours, sigma,
            num_viscous_iterations, num_elastic_iterations);
      },
      "floating_features"_a, "corresponding_features"_a, "floating_faces"_a,
      "floating_flags"_a, "weights"_a, "num_smoothing_neighbours"_a = 10,
      "sigma"_a = 3.0f, "num_viscous_iterations"_a = 50,
      "num_elastic_iterations"_a = 50,
      nb::call_guard<nb::gil_scoped_release>());

  m.def(
      "downsample_mesh",
      [](Eigen::Ref<const meshmonk::FeatureMat> features,
         Eigen::Ref<const meshmonk::FacesMat> faces,
         Eigen::Ref<const meshmonk::VecDynFloat> flags, float downsample_ratio)
          -> std::tuple<meshmonk::FeatureMat, meshmonk::FacesMat,
                        meshmonk::VecDynFloat, Eigen::VectorXi> {
        return meshmonk::downsample_mesh(features, faces, flags,
                                         downsample_ratio);
      },
      "features"_a, "faces"_a, "flags"_a, "downsample_ratio"_a = 0.8f,
      nb::call_guard<nb::gil_scoped_release>());

  m.def(
      "scale_shift_mesh",
      [](Eigen::Ref<const meshmonk::FeatureMat> previous_features,
         Eigen::Ref<const Eigen::VectorXi> previous_indices,
         Eigen::Ref<const Eigen::VectorXi> new_indices)
          -> meshmonk::FeatureMat {
        return meshmonk::scale_shift_mesh(previous_features, previous_indices,
                                          new_indices);
      },
      "previous_features"_a, "previous_indices"_a, "new_indices"_a,
      nb::call_guard<nb::gil_scoped_release>());

  m.def(
      "compute_normals",
      [](Eigen::Ref<const meshmonk::Vec3Mat> positions,
         Eigen::Ref<const meshmonk::FacesMat> faces) -> meshmonk::Vec3Mat {
        return meshmonk::compute_normals(positions, faces);
      },
      "positions"_a, "faces"_a, nb::call_guard<nb::gil_scoped_release>());

  m.def(
      "set_log_level",
      [](meshmonk::LogLevel level) { meshmonk::set_log_level(level); },
      "level"_a);

  // -----------------------------------------------------------------------
  // Profiling bindings (always compiled — no-ops when MESHMONK_PROFILING OFF)
  // -----------------------------------------------------------------------
  m.def("profiling_reset", []() {
#ifdef MESHMONK_PROFILING
    g_profiler.reset();
#endif
  });

  m.def("profiling_peek",
    []() -> std::unordered_map<std::string,
                                std::unordered_map<std::string, uint64_t>> {
#ifdef MESHMONK_PROFILING
      auto raw = g_profiler.snapshot();
      std::unordered_map<std::string,
                         std::unordered_map<std::string, uint64_t>> out;
      for (auto& [k, v] : raw) {
        out[k] = {{"total_us", v.total_us}, {"count", v.count}};
      }
      return out;
#else
      return {};
#endif
    });

  m.def("profiling_dump",
    []() -> std::unordered_map<std::string,
                                std::unordered_map<std::string, uint64_t>> {
#ifdef MESHMONK_PROFILING
      auto raw = g_profiler.snapshot();
      std::unordered_map<std::string,
                         std::unordered_map<std::string, uint64_t>> out;
      for (auto& [k, v] : raw) {
        out[k] = {{"total_us", v.total_us}, {"count", v.count}};
      }
      g_profiler.reset();
      return out;
#else
      return {};
#endif
    });

  m.def("profiling_enabled", []() -> bool {
#ifdef MESHMONK_PROFILING
    return true;
#else
    return false;
#endif
  });

  m.def("profiling_calibrate", [](size_t n) -> uint64_t {
#ifdef MESHMONK_PROFILING
    return g_profiler.calibrate(n);
#else
    return 0;
#endif
  }, nb::arg("n") = static_cast<size_t>(1'000'000));
}
