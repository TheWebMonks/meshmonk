#include "CorrespondenceFilter.hpp"


namespace registration {


void CorrespondenceFilter::set_floating_input(const FeatureMat * const inFloatingFeatures,
                                              const VecDynFloat * const inFloatingFlags)
{
    //# Set input
    _inFloatingFeatures = inFloatingFeatures;
    _inFloatingFlags = inFloatingFlags;

    //# Update internal parameters
    _numFloatingElements = _inFloatingFeatures->rows();
    _numAffinityElements = _numFloatingElements * _numNeighbours;

    //# Update the neighbour finder
    _neighbourFinder.set_queried_points(_inFloatingFeatures);
}

void CorrespondenceFilter::set_target_input(const FeatureMat * const inTargetFeatures,
                                            const VecDynFloat * const inTargetFlags)
{
    //# Set input
    _inTargetFeatures = inTargetFeatures;
    _inTargetFlags = inTargetFlags;

    //# Update internal parameters
    _numTargetElements = _inTargetFeatures->rows();
    _numAffinityElements = _numFloatingElements * _numNeighbours;

    //# Update the neighbour finder
    _neighbourFinder.set_source_points(_inTargetFeatures);
}


void CorrespondenceFilter::set_parameters(const size_t numNeighbours,
                                          const float flagThreshold)
{
    _numNeighbours = numNeighbours;
    _flagThreshold = flagThreshold;
    _numAffinityElements = _numFloatingElements * _numNeighbours;
    _neighbourFinder.set_parameters(_numNeighbours);
}


void CorrespondenceFilter::set_affinity_normalization(const bool normalizeAffinity){
    _normalizeAffinity = normalizeAffinity;
}

void CorrespondenceFilter::_update_affinity() {
    /*
    # GOALthe
    For each element in _inFloatingFeatures, we're going to determine affinity weights
    which link it to elements in _inTargetFeatures. This is based on the Euclidean
    distance of k nearest neighbours found for each element of _inFloatingFeatures.
    */

    //# Initialization
    //## Initialize the sparse affinity matrix
    _affinity = SparseMat(_numFloatingElements, _numTargetElements);
    _affinity.reserve(_numAffinityElements);
    //## Initialize a vector of triplets which is used to insert elements into
    //## the affinity matrix later.
    std::vector<Triplet> affinityElements(_numAffinityElements, Triplet(0,0,0.0f));

    //## Obtain pointers to neighbouring indices and (squared) distances:
    const MatDynInt neighbourIndices = _neighbourFinder.get_indices();
    const MatDynFloat neighbourSquaredDistances = _neighbourFinder.get_distances();

    //# Compute the affinity matrix
    //## Loop over the first feature set to determine their affinity with the
    //## second set.
    size_t i = 0;
    size_t j = 0;
    Vec3Float floatingNormal = Vec3Float::Zero();
    Vec3Float targetNormal = Vec3Float::Zero();
    unsigned int counter = 0;
    for ( ; i < _numFloatingElements ; i++) {
        floatingNormal = _inFloatingFeatures->row(i).tail(3);
        //### Loop over each found neighbour
        for ( j = 0 ; j < _numNeighbours ; j++) {
            //### Get index of neighbour and squared distance to it
            const int neighbourIndex = neighbourIndices(i,j);
            float distanceSquared = neighbourSquaredDistances(i,j);

            //### For numerical stability, check if the distance is very small
            const float eps1 = 0.000001f;
            if (distanceSquared < eps1) {distanceSquared = eps1;}

            //### Compute the affinity element as 1/distance*distance
            float affinityElement = 1.0f / distanceSquared;

            //### Incorporate the orientation
            targetNormal = _inTargetFeatures->row(neighbourIndex).tail(3);
            float dotProduct = floatingNormal.dot(targetNormal);
            float orientationWeight = dotProduct / 2.0f + 0.5f;
            affinityElement *= orientationWeight;

            //### Check for numerical stability (the affinity elements will be
            //### normalized later, so dividing by a sum of tiny elements might
            //### go wrong.
            const float eps2 = 0.0001f;
            if (affinityElement < eps2) {affinityElement = eps2;}

            //### Write result to the affinity triplet list
            affinityElements[counter] = Triplet(i,neighbourIndex,affinityElement);
            counter++;
        }
    }
    //## Construct the sparse matrix with the computed element list
    _affinity.setFromTriplets(affinityElements.begin(),
                                affinityElements.end());

    //# Normalize the rows of the affinity matrix
    if (_normalizeAffinity) {
        normalize_sparse_matrix(_affinity);
    }

}//end wknn_affinity()


void CorrespondenceFilter::update() {

    //# Update the neighbour indices and distances
    _neighbourFinder.set_queried_points(_inFloatingFeatures);
    _neighbourFinder.update();

    //# Update the (sparse) affinity matrix
    _update_affinity();

    if (_ioCorrespondingFeatures != NULL) {
        //# Use the affinity weights to determine corresponding features and flags.
        if (_normalizeAffinity) {
            BaseCorrespondenceFilter::_affinity_to_correspondences();
        }
        else {
            std::cout << "Can't compute correspondences without normalizing the affinity matrix!" << std::endl;
        }
    }
    else {
        /*
        affinity is computed and get be requested from the filter, but output correspondences
        are not generated if the filter doesn't know where to write the output.
        */
    }
}//end wkkn_correspondences()

}
