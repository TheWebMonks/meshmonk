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


void SymmetricCorrespondenceFilter::update() {

    //# Update the I/O for the push- and pull-filters
    _pushFilter.set_floating_input(_inFloatingFeatures, _inFloatingFlags);
    _pushFilter.set_target_input(_inTargetFeatures, _inTargetFlags);
    _pullFilter.set_floating_input(_inTargetFeatures, _inTargetFlags);
    _pullFilter.set_target_input(_inFloatingFeatures, _inFloatingFlags);


    //# Update the (sparse) affinity matrix
    _update_affinity();

    //# Use the affinity weights to determine corresponding features and flags.
    BaseCorrespondenceFilter::_affinity_to_correspondences();

}//end wkkn_correspondences()

}//namespace registration
