#ifndef NEIGHBOURFINDER_HPP
#define NEIGHBOURFINDER_HPP

#include <Eigen/Dense>
#include <nanoflann.hpp>

typedef Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic> MatDynInt; //matrix MxN of type unsigned int
typedef Eigen::Matrix< float, Eigen::Dynamic, Eigen::Dynamic> MatDynFloat;

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
        NeighbourFinder();
        ~NeighbourFinder(); //destructor

        void set_source_points(const VecMatType * const inQueriedPoints);
        void set_queried_points(const VecMatType * const _inSourcePoints);
        const MatDynInt * const get_indices() { return _outNeighbourIndices;}
        const MatDynFloat * const get_distances() { return _outNeighbourSquaredDistances;}
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
        nanoflann::KDTreeEigenMatrixAdaptor<VecMatType> _kdTree;

        //# Interal parameters
        size_t _numDimensions = 0;
        size_t _numSourceElements = 0;
        size_t _numQueriedElements = 0;
        size_t _numNeighbours = 3;
        size_t _leafSize = 15;

};


} //namespace registration

#endif // NEIGHBOURFINDER_HPP
