#include "RigidRegistration.hpp"
#include "meshmonk/profiling.hpp"
#include <memory>

namespace registration {

void RigidRegistration::set_input(FeatureMat *const ioFloatingFeatures,
                                  const FeatureMat *const inTargetFeatures,
                                  const VecDynFloat *const inFloatingFlags,
                                  const VecDynFloat *const inTargetFlags) {
  _ioFloatingFeatures = ioFloatingFeatures;
  _inTargetFeatures = inTargetFeatures;
  _inFloatingFlags = inFloatingFlags;
  _inTargetFlags = inTargetFlags;
} // end set_input()

void RigidRegistration::set_parameters(bool symmetric, size_t numNeighbours,
                                       float flagThreshold,
                                       bool equalizePushPull, float kappaa,
                                       bool inlierUseOrientation,
                                       size_t numIterations, bool useScaling) {
  _symmetric = symmetric;
  _numNeighbours = numNeighbours;
  _flagThreshold = flagThreshold;
  _equalizePushPull = equalizePushPull;
  _kappaa = kappaa;
  _inlierUseOrientation = inlierUseOrientation;
  _numIterations = numIterations;
  _useScaling = useScaling;
} // end set_parameters()

void RigidRegistration::update() {
#ifdef MESHMONK_PROFILING
  auto _t = g_profiler.scoped("RigidRegistration::update");
#endif

  // # Initializes
  size_t numFloatingVertices = _ioFloatingFeatures->rows();
  FeatureMat correspondingFeatures =
      FeatureMat::Zero(numFloatingVertices, registration::NUM_FEATURES);
  VecDynFloat correspondingFlags = VecDynFloat::Zero(numFloatingVertices);

  // # Set up the filters
  // ## Correspondence Filter
  std::unique_ptr<BaseCorrespondenceFilter> correspondenceFilter;

  if (_symmetric) {
    auto *f = new SymmetricCorrespondenceFilter();
    f->set_parameters(_numNeighbours, _flagThreshold, _equalizePushPull);
    correspondenceFilter.reset(f);
  } else {
    auto *f = new CorrespondenceFilter();
    f->set_parameters(_numNeighbours, _flagThreshold);
    correspondenceFilter.reset(f);
  }
  correspondenceFilter->set_floating_input(_ioFloatingFeatures,
                                           _inFloatingFlags);
  correspondenceFilter->set_target_input(_inTargetFeatures, _inTargetFlags);
  correspondenceFilter->set_output(&correspondingFeatures, &correspondingFlags);

  // ## Inlier Filter
  VecDynFloat floatingWeights = VecDynFloat::Ones(numFloatingVertices);
  InlierDetector inlierDetector;
  inlierDetector.set_input(_ioFloatingFeatures, &correspondingFeatures,
                           &correspondingFlags);
  inlierDetector.set_output(&floatingWeights);
  inlierDetector.set_parameters(_kappaa, _inlierUseOrientation);
  // ## Transformation Filter
  RigidTransformer rigidTransformer;
  rigidTransformer.set_input(&correspondingFeatures, &floatingWeights);
  rigidTransformer.set_output(_ioFloatingFeatures);
  rigidTransformer.set_parameters(_useScaling);

  // # Perform ICP
  for (size_t iteration = 0; iteration < _numIterations; iteration++) {
    // # Correspondences
    correspondenceFilter->set_floating_input(_ioFloatingFeatures,
                                             _inFloatingFlags);
    correspondenceFilter->set_target_input(_inTargetFeatures, _inTargetFlags);
    correspondenceFilter->update();

    // # Inlier Detection
    inlierDetector.update();

    // # Transformation
    rigidTransformer.update();

    // # Update final transformation matrix
    Mat4Float currentTransform = rigidTransformer.get_transformation();
    _transformationMatrix = currentTransform * _transformationMatrix;
  }

  // correspondenceFilter is automatically deleted by unique_ptr

} // end update()

} // namespace registration
