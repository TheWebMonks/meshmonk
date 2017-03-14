#ifndef SYMMETRICCORRESPONDENCEFILTER_HPP
#define SYMMETRICCORRESPONDENCEFILTER_HPP

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <stdio.h>
#include "../global.hpp"
#include "BaseCorrespondenceFilter.hpp"
#include "CorrespondenceFilter.hpp"
#include "helper_functions.hpp"

typedef Eigen::VectorXf VecDynFloat;
typedef Eigen::Matrix< float, Eigen::Dynamic, registration::NUM_FEATURES> FeatureMat; //matrix Mx6 of type float
typedef Eigen::Matrix< float, 1, registration::NUM_FEATURES> FeatureVec; //matrix Mx6 of type float
typedef Eigen::Vector3f Vec3Float;
typedef Eigen::SparseMatrix<float, 0, int> SparseMat;
typedef Eigen::Triplet<float> Triplet;
typedef Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic> MatDynInt; //matrix MxN of type unsigned int
typedef Eigen::Matrix< float, Eigen::Dynamic, Eigen::Dynamic> MatDynFloat;

namespace registration {

class SymmetricCorrespondenceFilter: public BaseCorrespondenceFilter
{
    /*
    # GOAL
    For each element in inFloatingFeatures, we're going to find corresponding
    features in the elements of inTargetFeatures.

    A push-pull approach is used to find correspondences. Not only
    will we look for correspondences for inFloatingFeatures in the
    inTargetFeatures set, we'll also do the opposite: find correspondences for
    the inTargetFeatures set in the inFloatingFeatures set. These findings are
    combined which creates an effect where the inTargetFeatures sort of attract
    the inFloatingFeatures towards itself.

    # INPUTS
    -inFloatingFeatures
    -inTargetFeatures
    -inTargetFlags

    # PARAMETERS
    -numNeighbours(=3):
    number of nearest neighbours

    # OUTPUT
    -outCorrespondingFeatures
    -outCorrespondingFlags
    */

    public:
        //CorrespondenceFilter(); //default constructor
        //~CorrespondenceFilter(); //destructor

        void set_floating_input(const FeatureMat * const inFloatingFeatures,
                                const VecDynFloat * const inFloatingFlags);
        void set_target_input(const FeatureMat * const inTargetFeatures,
                            const VecDynFloat * const inTargetFlags);
        void set_parameters(const size_t numNeighbours);
        void update();

    protected:

    private:

        //# Internal Data structures
        CorrespondenceFilter _pushFilter;
        CorrespondenceFilter _pullFilter;

        //# Internal functions
        //## Function to update the sparse affinity matrix
        void _update_affinity();
        //## Function to convert the sparse affinity weights into corresponding
        //## features and flags
        void _affinity_to_correspondences();
};

}//namespace registration

#endif // SYMMETRICCORRESPONDENCEFILTER_HPP
