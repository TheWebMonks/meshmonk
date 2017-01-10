#ifndef CORRESPONDENCEFILTER_HPP
#define CORRESPONDENCEFILTER_HPP

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <stdio.h>
#include "../global.hpp"
#include <NeighbourFinder.hpp>
#include <helper_functions.hpp>

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
        NeighbourFinder _neighbourFinder;
        SparseMat _affinity;


        //# Internal Parameters
        size_t _numFloatingElements = 0;
        size_t _numTargetElements = 0;
        size_t _numAffinityElements = 0;
        float _flagRoundingLimit = 0.9;

        //# Internal functions
        //## Function to update the sparse affinity matrix
        void _update_affinity();
        //## Function to convert the sparse affinity weights into corresponding
        //## features and flags
        void _affinity_to_correspondences();
};




}//namespace registration

#endif // CORRESPONDENCEFILTER_HPP
