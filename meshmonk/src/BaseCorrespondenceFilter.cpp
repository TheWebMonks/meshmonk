#include "BaseCorrespondenceFilter.hpp"

namespace registration {

BaseCorrespondenceFilter::BaseCorrespondenceFilter()
{
    //ctor
}

BaseCorrespondenceFilter::~BaseCorrespondenceFilter()
{
    //dtor
}


void BaseCorrespondenceFilter::set_output(FeatureMat * const ioCorrespondingFeatures,
                                    VecDynFloat * const ioCorrespondingFlags)
{
    _ioCorrespondingFeatures = ioCorrespondingFeatures;
    _ioCorrespondingFlags = ioCorrespondingFlags;

}//end set_output

void BaseCorrespondenceFilter::_affinity_to_correspondences(){
    /*
    # GOAL
    This function computes all the corresponding features and flags,
    when the affinity for a set of features and flags is given.

    # PARAMETERS
    -_flagRoundingLimit:
    Flags are binary. Anything over this double is rounded up, whereas anything
    under it is rounded down. A suggested cut-off is 0.9, which means that if
    an element flagged zero contributes 10 percent or to the affinity, the
    corresponding element should be flagged a zero as well.

    # RESULT
    -outCorrespondingFeatures
    -outCorrespondingFlags
    */

    //# Simple computation of corresponding features and flags
    *_ioCorrespondingFeatures = _affinity * (*_inTargetFeatures);
    *_ioCorrespondingFlags = _affinity * (*_inTargetFlags);

    //# Flag correction.
    //## Flags are binary. We will round them down if lower than the flag
    //## rounding limit (see explanation in parameter description).
    if (_flagThreshold >= 1.0f) {std::cerr << "corresponding flag threshold equals " << _flagThreshold << " but has to be between 0.0 and 1.0!" <<std::endl;}
    for (size_t i = 0 ; i < _numFloatingElements ; i++) {
        if ((*_ioCorrespondingFlags)[i] > _flagThreshold){
            (*_ioCorrespondingFlags)[i] = 1.0;
        }
        else {
            (*_ioCorrespondingFlags)[i] = 0.0;
        }
    }
}

}//namespace registration
