#ifndef INLIERDETECTOR_HPP
#define INLIERDETECTOR_HPP

#include <Eigen/Dense>
#include <stdio.h>
#include <iostream>
#include <map>
#include "../global.hpp"
#include "NeighbourFinder.hpp"

typedef Eigen::VectorXf VecDynFloat;
typedef Eigen::Matrix< float, Eigen::Dynamic, registration::NUM_FEATURES> FeatureMat; //matrix Mx6 of type float
typedef Eigen::Matrix< float, 1, registration::NUM_FEATURES> FeatureVec; //matrix Mx6 of type float
typedef Eigen::Vector3f Vec3Float;
typedef Eigen::Matrix< float, Eigen::Dynamic, 3> Vec3Mat; //matrix Mx3 of type float
typedef Eigen::Matrix< float, Eigen::Dynamic, Eigen::Dynamic> MatDynFloat;

namespace registration {

class InlierDetector
{
    /*
    # GOAL
    Determine which elements are inliers ('normal') and which are outliers
    ('abnormal'). There are different ways to to that, but they should all
    result in a scalar assigmnent that represents the probability of an element
    being an inlier.

    # INPUTS
    -_inFeatures
    -_inCorrespondingFeatures
    -_inCorrespondingFlags

    # PARAMETERS
    -_kappa(=3): Mahalanobis distance that determines cut-off in- vs outliers

    # OUTPUTS
    -_ioProbability

    # RETURNS
    */

    private:
        //# Inputs
        const FeatureMat * _inFeatures = NULL;
        const FeatureMat * _inCorrespondingFeatures = NULL;
        const VecDynFloat * _inCorrespondingFlags = NULL;

        //# Outputs
        VecDynFloat *_ioProbability = NULL;

        //# User Parameters
        float _kappa = 3.0;
        bool _useOrientation = true;

        //# Internal variables
        size_t _numElements;
        float _minimalSigma = 0.1f;
        float _maximalSigma = 10.0f;
        const size_t _numNeighbours = 10;
        const float _minWeight = 0.00001f;
        const size_t _numSmoothingPasses = 2;

        //# Internal Data structures
        NeighbourFinder<Vec3Mat> _neighbourFinder;
        MatDynFloat _smoothingWeights;

        //# Internal functions
        //## Find nearest neighbours (required for smoothing inlier weights)
        void _determine_neighbours();
        //## Update the weights used for smoothing the inlier weights (distance-based)
        void _update_smoothing_weights();
        //## Smooth the inlier weights
        void _smooth_inlier_weights();

    protected:

    public:
        void set_input(const FeatureMat * const inFeatures,
                        const FeatureMat * const inCorrespondingFeatures,
                        const VecDynFloat * const inCorrespondingFlags);
        void set_output(VecDynFloat * const _ioProbability);
        void set_parameters(const float kappa, const bool useOrientation);
        void update();
};

} //namespace registration

#endif
