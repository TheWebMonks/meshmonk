#include "CorrespondenceFilter.hpp"


namespace registration {

CorrespondenceFilter::CorrespondenceFilter()
{
    //ctor
}

CorrespondenceFilter::~CorrespondenceFilter()
{
    //destructor
    if (affinity != NULL) { delete affinity; affinity = NULL;}
}


void CorrespondenceFilter::set_input(const FeatureMat * const inFloatingFeatures,
                                    const FeatureMat * const inTargetFeatures,
                                    const VecDynFloat * const inTargetFlags)
{
    //# Set input
    _inFloatingFeatures = inFloatingFeatures;
    _inTargetFeatures = inTargetFeatures;
    _inTargetFlags = inTargetFlags;

    //# Update internal parameters
    _numFloatingElements = _inFloatingFeatures->rows();
    _numTargetElements = _inTargetFeatures->rows();
    _numAffinityElements = _numFloatingElements * _numNeighbours;

    //# Update internal data structures
    if (affinity != NULL) { delete affinity; affinity = NULL;}
    affinity = new SparseMat(_numFloatingElements, _numTargetElements);
    affinity->reserve(_numAffinityElements);
}

void CorrespondenceFilter::set_output(FeatureMat * const ioCorrespondingFeatures,
                                    VecDynFloat * const ioCorrespondingFlags)
{
    _ioCorrespondingFeatures = ioCorrespondingFeatures;
    _ioCorrespondingFlags = ioCorrespondingFlags;
}

void CorrespondenceFilter::set_parameters(const size_t numNeighbours)
{
    _numNeighbours = numNeighbours;
}

void CorrespondenceFilter::update() {

//    //# Compute the affinity for inFloatingFeatures towards inTargetFeatures
//    wknn_affinity(inFloatingFeatures, inTargetFeatures, affinity, paramK);
//
//    //# Ccompute symmetric correspondences (if required)
//    if (paramSymmetric == true) {
//        affinityPull = SparseMat(numTargetVertices, numFloatingVertices);
//        //## Compute the affinity for inTargetFeatures towards inFloatingFeatures
//        wknn_affinity(inTargetFeatures, inFloatingFeatures, affinityPull, paramK);
//
//        //## Fuse the two affinities together.
//        fuse_affinities(affinity, affinityPull);
//    }
//
//    //## Compute corresponding features as affinity matrix multiplied with the
//    //## target features.
//    affinity_to_correspondences(inTargetFeatures, inTargetFlags, affinity,
//                                outCorrespondingFeatures, outCorrespondingFlags,
//                                0.9);
}//end wkkn_correspondences()

}
