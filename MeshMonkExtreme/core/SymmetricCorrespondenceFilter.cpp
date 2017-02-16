#include "SymmetricCorrespondenceFilter.hpp"

namespace registration {


void SymmetricCorrespondenceFilter::set_floating_input(const FeatureMat * const inFloatingFeatures,
                                                       const VecDynFloat * const inFloatingFlags)
{
    //# Set input
    _inFloatingFeatures = inFloatingFeatures;
    _inFloatingFlags = inFloatingFlags;

    //# Update internal parameters
    _numFloatingElements = _inFloatingFeatures->rows();

    //# Update the push and pull filters
    _pushFilter.set_floating_input(_inFloatingFeatures);
    _pullFilter.set_target_input(_inFloatingFeatures, _inFloatingFlags);
}

void SymmetricCorrespondenceFilter::set_target_input(const FeatureMat * const inTargetFeatures,
                                            const VecDynFloat * const inTargetFlags)
{
    //# Set input
    _inTargetFeatures = inTargetFeatures;
    _inTargetFlags = inTargetFlags;

    //# Update internal parameters
    _numTargetElements = _inTargetFeatures->rows();

    //# Update the push and pull filters
    _pushFilter.set_target_input(_inTargetFeatures, _inTargetFlags);
    _pullFilter.set_floating_input(_inTargetFeatures);
}

void SymmetricCorrespondenceFilter::set_output(FeatureMat * const ioCorrespondingFeatures,
                                    VecDynFloat * const ioCorrespondingFlags)
{
    _ioCorrespondingFeatures = ioCorrespondingFeatures;
    _ioCorrespondingFlags = ioCorrespondingFlags;
}

void SymmetricCorrespondenceFilter::set_parameters(const size_t numNeighbours)
{
    _numNeighbours = numNeighbours;
    _pushFilter.set_parameters(_numNeighbours);
    _pullFilter.set_parameters(_numNeighbours);
}

void SymmetricCorrespondenceFilter::_update_affinity() {

    //# Compute the push and pull affinity
    _pushFilter.update();
    _pullFilter.update();

    //# Get the affinities
    _affinity = _pushFilter.get_affinity();
    const SparseMat pullAffinity = _pullFilter.get_affinity();

    //# Fuse the affinities
    fuse_affinities(_affinity, pullAffinity); //helper function to combine affinity matrices

}//end wknn_affinity()

void SymmetricCorrespondenceFilter::_affinity_to_correspondences(){
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
    for (size_t i = 0 ; i < _numFloatingElements ; i++) {
        if ((*_ioCorrespondingFlags)[i] > _flagRoundingLimit){
            (*_ioCorrespondingFlags)[i] = 1.0;
        }
        else {
            (*_ioCorrespondingFlags)[i] = 0.0;
        }
    }
}

void SymmetricCorrespondenceFilter::update() {

    //# Update the I/O for the push- and pull-filters
    _pushFilter.set_floating_input(_inFloatingFeatures);
    _pushFilter.set_target_input(_inTargetFeatures, _inTargetFlags);
    _pullFilter.set_floating_input(_inTargetFeatures);
    _pullFilter.set_target_input(_inFloatingFeatures, _inFloatingFlags);


    //# Update the (sparse) affinity matrix
    _update_affinity();

    //# Use the affinity weights to determine corresponding features and flags.
    _affinity_to_correspondences();

}//end wkkn_correspondences()

}//namespace registration
