#include "CorrespondenceFilter.hpp"


namespace registration {


void CorrespondenceFilter::set_floating_input(const FeatureMat * const inFloatingFeatures)
{
    //# Set input
    _inFloatingFeatures = inFloatingFeatures;

    //# Update internal parameters
    _numFloatingElements = _inFloatingFeatures->rows();
    _numAffinityElements = _numFloatingElements * _numNeighbours;

    //# Update the neighbour finder
    //NeighbourFinder<FeatureMat> testje = NeighbourFinder<FeatureMat>();
    //testje.set_queried_points(_inFloatingFeatures);
    //_neighbourFinder.set_queried_points(_inFloatingFeatures);
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
    //_neighbourFinder.set_source_points(_inTargetFeatures);
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

void CorrespondenceFilter::_update_affinity() {
    /*
    # GOAL
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
    unsigned int counter = 0;
    for ( ; i < _numFloatingElements ; i++) {
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
            const float eps2 = 0.0001f;
            //### Check for numerical stability (the affinity elements will be
            //### normalized later, so dividing by a sum of tiny elements might
            //### go wrong.
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
    normalize_sparse_matrix(_affinity);

}//end wknn_affinity()

void CorrespondenceFilter::_affinity_to_correspondences(){
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

void CorrespondenceFilter::update() {

    //# Update the neighbour indices and distances
    _neighbourFinder.update();

    //# Update the (sparse) affinity matrix
    _update_affinity();

    //# Use the affinity weights to determine corresponding features and flags.
    _affinity_to_correspondences();


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
