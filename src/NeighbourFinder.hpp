#ifndef NEIGHBOURFINDER_HPP
#define NEIGHBOURFINDER_HPP

#include <Eigen/Dense>
#include <nanoflann.hpp>
#include "../global.hpp"

typedef Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic> MatDynInt; //matrix MxN of type unsigned int
typedef Eigen::Matrix< float, Eigen::Dynamic, Eigen::Dynamic> MatDynFloat;
typedef Eigen::Matrix< float, Eigen::Dynamic, registration::NUM_FEATURES> FeatureMat;
typedef Eigen::Matrix< float, Eigen::Dynamic, 3> Vec3Mat;

namespace registration {

template <typename VecMatType>
class NeighbourFinder
{
    /*
    GOAL
    This class searches for the k nearest neighbours in 'inSourcePoints' for
    each element in the 'inQueriedPoints' set. It outputs the indices of each
    neighbour and the squared (!) distances between each element of
    'inQueriedPoints' and its neighbours.

    INPUT
    -inQueriedPoints:
    -inSourcePoints:

    PARAMETERS
    -numNeighbours(= 3): number of nearest neighbours
    -leafSize(= 15): should be between 5 and 50 or so

    OUTPUT
    -outNeighbourIndices
    -outNeighbourSquaredDistances.
    */

    public:
        //NeighbourFinder();
        ~NeighbourFinder(); //destructor

        void set_source_points(const VecMatType * const inQueriedPoints);
        void set_queried_points(const VecMatType * const _inSourcePoints);
        MatDynInt get_indices() const { return _outNeighbourIndices;}
        MatDynFloat get_distances() const { return _outNeighbourSquaredDistances;}
        void set_parameters(const size_t numNeighbours);
        void update();

    protected:

    private:
        //# Inputs
        const VecMatType * _inQueriedPoints = NULL;
        const VecMatType * _inSourcePoints = NULL;

        //# Outputs
        MatDynInt _outNeighbourIndices;
        MatDynFloat _outNeighbourSquaredDistances;
        //# User parameters

        //# Internal Data structures
        nanoflann::KDTreeEigenMatrixAdaptor<VecMatType> * _kdTree = NULL;

        //# Interal parameters
        size_t _numDimensions = 0;
        size_t _numSourceElements = 0;
        size_t _numQueriedElements = 0;
        size_t _numNeighbours = 3;
        size_t _leafSize = 15;

};

/*
See http://stackoverflow.com/questions/495021/why-can-templates-only-be-implemented-in-the-header-file
for an explanation why we instantiate our templated class with every matrix type we will use here.
*/
//
template class NeighbourFinder<FeatureMat>;
template class NeighbourFinder<Vec3Mat>;
//


template <typename VecMatType>
NeighbourFinder<VecMatType>::~NeighbourFinder()
{
    //destructor
    if (_kdTree != NULL) { delete _kdTree; _kdTree = NULL;}
}

template <typename VecMatType>
void NeighbourFinder<VecMatType>::set_source_points(const VecMatType * const inSourcePoints){
    //# Set input
    _inSourcePoints = inSourcePoints;

    //# Update internal parameters
    _numDimensions = _inSourcePoints->cols();
    _numSourceElements = _inSourcePoints->rows();

    //# Update internal data structures
    //## The kd-tree has to be rebuilt.
    if (_kdTree != NULL) { delete _kdTree; _kdTree = NULL;}
    _kdTree = new nanoflann::KDTreeEigenMatrixAdaptor<VecMatType>(_numDimensions,
                                                                *_inSourcePoints,
                                                                _leafSize);
    _kdTree->index->buildIndex();
}


template <typename VecMatType>
void NeighbourFinder<VecMatType>::set_queried_points(const VecMatType * const inQueriedPoints){
    //# Set input
    _inQueriedPoints = inQueriedPoints;

    //# Update internal parameters
    _numQueriedElements = _inQueriedPoints->rows();

    //# Adjust internal data structures
    //## The indices and distance matrices have to be resized.
    _outNeighbourIndices.setZero(_numQueriedElements,_numNeighbours);
    _outNeighbourSquaredDistances.setZero(_numQueriedElements,_numNeighbours);

}


template <typename VecMatType>
void NeighbourFinder<VecMatType>::set_parameters(const size_t numNeighbours){
    //# Check if what user requests, changes the parameter value
    bool parameterChanged = false;
    if (_numNeighbours != numNeighbours){ parameterChanged = true;}

    //# Set parameter
    _numNeighbours = numNeighbours;

    //# Resize the output matrices if the parameter is changed
    if (parameterChanged == true) {
        _outNeighbourIndices.setZero(_numQueriedElements,_numNeighbours);
        _outNeighbourSquaredDistances.setZero(_numQueriedElements,_numNeighbours);
    }
}

template <typename VecMatType>
void NeighbourFinder<VecMatType>::update(){

    //# Query the kd-tree
    //## Loop over the queried features
    //### Initialize variables we'll need during the loop
    unsigned int i = 0;
    unsigned int j = 0;
    std::vector<float> queriedFeature(_numDimensions);
    std::vector<size_t> neighbourIndices(_numNeighbours);
    std::vector<float> neighbourSquaredDistances(_numNeighbours);
    nanoflann::KNNResultSet<float> knnResultSet(_numNeighbours);

    //### Execute loop
    for ( ; i < _numQueriedElements ; ++i ) {
        //### Initiliaze the knnResultSet
        knnResultSet.init(&neighbourIndices[0], &neighbourSquaredDistances[0]);

        //### convert input features to 'queriedFeature' std::vector structure
        //### (required by nanoflann's kd-tree).
        for (j = 0 ; j < _numDimensions ; ++j) {
            queriedFeature[j] = (*_inQueriedPoints)(i,j);
        }

        //### Query the kd-tree
//        size_t numNeighboursFound = kdTree.knnSearch(&queriedFeature[0], _numNeighbours, &neighbourIndices[0], &neighbourSquaredDistances[0]);
        _kdTree->index->findNeighbors(knnResultSet, &queriedFeature[0],
                                    nanoflann::SearchParams(32, 0.0001 /*eps*/, true));

        //### Copy the result into the outputs by looping over the k nearest
        //### neighbours
        for (j = 0 ; j < _numNeighbours ; ++j) {
            _outNeighbourIndices(i,j) = neighbourIndices[j];
            _outNeighbourSquaredDistances(i,j) = neighbourSquaredDistances[j];
        }
    }
}//end k_nearest_neighbours()

} //namespace registration

#endif // NEIGHBOURFINDER_HPP
