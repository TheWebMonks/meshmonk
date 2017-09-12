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
    _pushFilter.set_floating_input(_inFloatingFeatures, _inFloatingFlags);
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
    _pullFilter.set_floating_input(_inTargetFeatures, _inTargetFlags);
}

void SymmetricCorrespondenceFilter::set_parameters(const size_t numNeighbours,
                                                   const float flagThreshold,
                                                   const bool equalizePushPull)
{
    _numNeighbours = numNeighbours;
    _flagThreshold = flagThreshold;
    _equalizePushPull = equalizePushPull;
    _pushFilter.set_parameters(_numNeighbours, _flagThreshold);
    _pullFilter.set_parameters(_numNeighbours, _flagThreshold);
    _pushFilter.set_affinity_normalization(false);
    _pullFilter.set_affinity_normalization(false);
}

void SymmetricCorrespondenceFilter::_update_push_and_pull() {

    //# Compute the push and pull affinity
    _pushFilter.update();
    _pullFilter.update();

    //# Get the affinities
    _affinity = _pushFilter.get_affinity();
    SparseMat pullAffinity = _pullFilter.get_affinity();

    //# Normalize the affinities before fusing them?
    if (_equalizePushPull) {
        normalize_sparse_matrix(_affinity);
        normalize_sparse_matrix(pullAffinity);
    }

    //# Fuse the affinities
    fuse_affinities(_affinity, pullAffinity); //helper function to combine affinity matrices

}//end wknn_affinity()


void SymmetricCorrespondenceFilter::update() {

    //# Update the I/O for the push- and pull-filters
    _pushFilter.set_floating_input(_inFloatingFeatures, _inFloatingFlags);
    _pushFilter.set_target_input(_inTargetFeatures, _inTargetFlags);
    _pullFilter.set_floating_input(_inTargetFeatures, _inTargetFlags);
    _pullFilter.set_target_input(_inFloatingFeatures, _inFloatingFlags);

    //# Set the output for the pull filter so that we can extract the pull flags later.
    FeatureMat pullFeatures = FeatureMat::Zero(_numTargetElements, NUM_FEATURES);
    VecDynFloat pullFlags = VecDynFloat::Ones(_numTargetElements);
    _pullFilter.set_output(&pullFeatures, &pullFlags);

    //# Update the push and pull filters.
    //## This gives us the affinity matrices for both push and pull directions,
    //## and the pull flags.
    _update_push_and_pull();

    //# The pull flags are updated. We now compute the corresponding pull flags.
    VecDynFloat correspondingPullFlags = VecDynFloat::Ones(_numFloatingElements);
    correspondingPullFlags = _affinity * pullFlags;

    //# Use the affinity weights to determine corresponding features and flags.
    BaseCorrespondenceFilter::_affinity_to_correspondences();

    //# Merge the just computed corresponding flags with the corresponding pull flags.
    (*BaseCorrespondenceFilter::_ioCorrespondingFlags) *= correspondingPullFlags;

}//end wkkn_correspondences()

}//namespace registration
