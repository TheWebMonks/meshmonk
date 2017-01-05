#include "NeighbourFinder.hpp"

namespace registration {

NeighbourFinder::NeighbourFinder()
{
    //ctor
}




void update(const VecMatType &inQueriedPoints,
                        const VecMatType &inSourcePoints,
                        MatDynInt &outNeighbourIndices,
                        MatDynFloat &outNeighbourSquaredDistances,
                        const size_t paramK = 3,
                        const size_t paramLeafsize = 15){
    /*
    GOAL
    This function searches for the k nearest neighbours in 'inSourcePoints' for
    each element in the 'inQueriedPoints' set. It outputs the indices of each
    neighbour and the squared (!) distances between each element of
    'inQueriedPoints' and its neighbours.

    INPUT
    -inQueriedPoints:
    -inSourcePoints:

    PARAMETERS
    -paramK(= 3): number of nearest neighbours
    -paramLeafsize(= 15): should be between 5 and 50 or so

    OUTPUT
    -outNeighbourIndices
    -outNeighbourSquaredDistances.
    */

    //# Info and Initialization
    const size_t dimension = inSourcePoints.cols();
    const size_t numSourceElements = inSourcePoints.rows();
    const size_t numQueriedElements = inQueriedPoints.rows();
    outNeighbourIndices = MatDynInt::Zero(numQueriedElements,paramK);
    outNeighbourSquaredDistances = MatDynFloat::Zero(numQueriedElements,paramK);

    //# Construct kd-tree
    nanoflann::KDTreeEigenMatrixAdaptor<VecMatType> kdTree(dimension, inSourcePoints, paramLeafsize);
    kdTree.index->buildIndex();

    //# Query the kd-tree
    //## Loop over the queried features
    //### Initialize variables we'll need during the loop
    unsigned int i = 0;
    unsigned int j = 0;
    std::vector<float> queriedFeature(dimension);
    std::vector<size_t> neighbourIndices(paramK);
    std::vector<float> neighbourSquaredDistances(paramK);
    nanoflann::KNNResultSet<float> knnResultSet(paramK);

    //### Execute loop
    for ( ; i < numQueriedElements ; ++i ) {
        //### Initiliaze the knnResultSet
        knnResultSet.init(&neighbourIndices[0], &neighbourSquaredDistances[0]);

        //### convert input features to 'queriedFeature' std::vector structure
        //### (required by nanoflann's kd-tree).
        for (j = 0 ; j < dimension ; ++j) {
            queriedFeature[j] = inQueriedPoints(i,j);
        }

        //### Query the kd-tree
//        size_t numNeighboursFound = kdTree.knnSearch(&queriedFeature[0], paramK, &neighbourIndices[0], &neighbourSquaredDistances[0]);
        kdTree.index->findNeighbors(knnResultSet, &queriedFeature[0],
                                    nanoflann::SearchParams(32, 0.0001 /*eps*/, true));

        //### Copy the result into the outputs by looping over the k nearest
        //### neighbours
        for (j = 0 ; j < paramK ; ++j) {
            outNeighbourIndices(i,j) = neighbourIndices[j];
            outNeighbourSquaredDistances(i,j) = neighbourSquaredDistances[j];
        }
    }
}//end k_nearest_neighbours()

}// namespace registration
