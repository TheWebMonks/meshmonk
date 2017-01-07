#ifndef CORRESPONDENCEFILTER_H
#define CORRESPONDENCEFILTER_H

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <stdio.h>
#include "../global.hpp"
#include <NeighbourFinder.hpp>

typedef Eigen::VectorXf VecDynFloat;
typedef Eigen::Matrix< float, Eigen::Dynamic, registration::NUM_FEATURES> FeatureMat; //matrix Mx6 of type float
typedef Eigen::Matrix< float, 1, registration::NUM_FEATURES> FeatureVec; //matrix Mx6 of type float
typedef Eigen::Vector3f Vec3Float;
typedef Eigen::SparseMatrix<float, 0, int> SparseMat;
typedef Eigen::Triplet<float> Triplet;
typedef Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic> MatDynInt; //matrix MxN of type unsigned int
typedef Eigen::Matrix< float, Eigen::Dynamic, Eigen::Dynamic> MatDynFloat;

namespace registration {

class CorrespondenceFilter
{

    /*
    # GOAL
    For each element in inFloatingFeatures, we're going to find corresponding
    features in the elements of inTargetFeatures.
    The used weight is one over the distance squared.

    # INPUTS
    -inFloatingFeatures
    -inTargetFeatures
    -inTargetFlags

    # PARAMETERS
    -numNeighbours(=3):
    number of nearest neighbours
    -paramSymmetric:
    If true, a push-pull approach is used to find correspondences. Not only
    will we look for correspondences for inFloatingFeatures in the
    inTargetFeatures set, we'll also do the opposite: find correspondences for
    the inTargetFeatures set in the inFloatingFeatures set. These findings are
    combined which creates an effect where the inTargetFeatures sort of attract
    the inFloatingFeatures towards itself.

    # OUTPUT
    -outCorrespondingFeatures
    -outCorrespondingFlags
    */

    public:
        CorrespondenceFilter(); //default constructor
        ~CorrespondenceFilter(); //destructor

        void set_floating_input(const FeatureMat * const inFloatingFeatures);
        void set_target_input(const FeatureMat * const inTargetFeatures,
                            const VecDynFloat * const inTargetFlags);
        void set_output(FeatureMat * const ioCorrespondingFeatures,
                        VecDynFloat * const ioCorrespondingFlags);
        void set_parameters(const size_t numNeighbours);
        void update();

    protected:

    private:
        //# Inputs
        const FeatureMat * _inFloatingFeatures = NULL;
        const FeatureMat * _inTargetFeatures = NULL;
        const VecDynFloat * _inTargetFlags = NULL;

        //# Outputs
        FeatureMat * _ioCorrespondingFeatures = NULL;
        VecDynFloat * _ioCorrespondingFlags = NULL;

        //# User Parameters
        size_t _numNeighbours = 3;

        //# Internal Data structures
        NeighbourFinder<FeatureMat> _neighbourFinder;
        SparseMat _affinity = NULL;


        //# Internal Parameters
        size_t _numFloatingElements = 0;
        size_t _numTargetElements = 0;
        size_t _numAffinityElements = 0;

        //# Internal functions
        void _update_affinity();
};




}//namespace registration

#endif // CORRESPONDENCEFILTER_H
