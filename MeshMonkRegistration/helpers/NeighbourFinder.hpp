#ifndef NEIGHBOURFINDER_HPP
#define NEIGHBOURFINDER_HPP

#include <Eigen/Dense>

typedef Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic> MatDynInt; //matrix MxN of type unsigned int
typedef Eigen::Matrix< float, Eigen::Dynamic, Eigen::Dynamic> MatDynFloat;

namespace registration {

template <typename VecMatType>
class NeighbourFinder
{
    public:
        NeighbourFinder();

        void set_input(const VecMatType * const _inQueriedPoints,
                        const VecMatType * const _inSourcePoints);
        const MatDynInt * const get_indices() { return _outNeighbourIndices;}
        const MatDynFloat * const get_distances() { return _outNeighbourSquaredDistances;}
        void set_parameters(const float kappa);
        void update();

    protected:

    private:
        //# Inputs
        const VecMatType * _inQueriedPoints = NULL;
        const VecMatType * _inSourcePoints = NULL;

        //# Outputs
        MatDynInt * _outNeighbourIndices = NULL;
        MatDynFloat * _outNeighbourSquaredDistances = NULL;
        //# User parameters

        //# Internal Data structures

        //# Interal parameters
        size_t _numNeighbours = 3;
        size_t _leafSize = 15;

};


} //namespace registration

#endif // NEIGHBOURFINDER_HPP
