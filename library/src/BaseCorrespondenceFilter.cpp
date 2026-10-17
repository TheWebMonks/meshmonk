#include "BaseCorrespondenceFilter.hpp"
#include "meshmonk/logger.hpp"

namespace registration {

BaseCorrespondenceFilter::BaseCorrespondenceFilter() {
  // ctor
}

BaseCorrespondenceFilter::~BaseCorrespondenceFilter() {
  // dtor
}

void BaseCorrespondenceFilter::set_output(
    FeatureMat *const ioCorrespondingFeatures,
    VecDynFloat *const ioCorrespondingFlags) {
  _ioCorrespondingFeatures = ioCorrespondingFeatures;
  _ioCorrespondingFlags = ioCorrespondingFlags;

} // end set_output

void BaseCorrespondenceFilter::_affinity_to_correspondences() {
  /*
  # GOAL
  This function computes all the corresponding features and flags,
  when the affinity for a set of features and flags is given.

  # PARAMETERS
  -_flagThreshold:
  Flags are binary. A weighted corresponding flag over this threshold
  (default 0.99, must be < 1) is rounded up to 1, anything else down to 0.

  # RESULT
  -outCorrespondingFeatures
  -outCorrespondingFlags
  */

  // # Simple computation of corresponding features and flags
  *_ioCorrespondingFeatures = _affinity * (*_inTargetFeatures);
  *_ioCorrespondingFlags = _affinity * (*_inTargetFlags);

  // # Flag correction.
  // ## Flags are binary. We will round them down if lower than the flag
  // ## rounding limit (see explanation in parameter description).
  if (_flagThreshold >= 1.0f) {
    meshmonk::log(meshmonk::LogLevel::Warning,
                  "corresponding flag threshold equals " +
                      std::to_string(_flagThreshold) +
                      " but has to be between 0.0 and 1.0!");
  }
  for (size_t i = 0; i < _numFloatingElements; i++) {
    if ((*_ioCorrespondingFlags)[i] > _flagThreshold) {
      (*_ioCorrespondingFlags)[i] = 1.0;
    } else {
      (*_ioCorrespondingFlags)[i] = 0.0;
    }
  }

  // # Merge corresponding and floating flags
  (*_ioCorrespondingFlags) =
      (*_ioCorrespondingFlags).cwiseProduct(*_inFloatingFlags);
}

} // namespace registration
