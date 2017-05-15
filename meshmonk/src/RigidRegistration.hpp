#ifndef RIGIDREGISTRATION_HPP
#define RIGIDREGISTRATION_HPP

#include <Eigen/Dense>
#include <stdio.h>
#include <memory.h>
#include <time.h>
#include "../global.hpp"
#include "CorrespondenceFilter.hpp"
#include "SymmetricCorrespondenceFilter.hpp"
#include "InlierDetector.hpp"
#include "RigidTransformer.hpp"

typedef Eigen::VectorXf VecDynFloat;
typedef Eigen::Matrix< float, Eigen::Dynamic, registration::NUM_FEATURES> FeatureMat; //matrix Mx6 of type float

namespace registration{

class RigidRegistration
{
    /*
    # GOAL
    This class performs icp-based rigid registration between two oriented pointclouds.

    # INPUTS
    -ioFloatingFeatures
    -inTargetFeatures
    -inFloatingFlags
    -inTargetFlags

    # PARAMETERS
    -numNeighbours(=3):
    number of nearest neighbours

    # OUTPUT
    -outCorrespondingFeatures
    -outCorrespondingFlags
    */

    public:

        void set_input(FeatureMat * const ioFloatingFeatures,
                       const FeatureMat * const inTargetFeatures,
                       const VecDynFloat * const inFloatingFlags,
                       const VecDynFloat * const inTargetFlags);
        void set_parameters(bool symmetric, size_t numNeighbours,
                            float flagThreshold, float kappaa,
                            bool inlierUseOrientation, size_t numIterations);

        void update();

    protected:

    private:
        //# Inputs/Outputs
        FeatureMat * _ioFloatingFeatures = NULL;
        const FeatureMat * _inTargetFeatures = NULL;
        const VecDynFloat * _inFloatingFlags = NULL;
        const VecDynFloat * _inTargetFlags = NULL;

        //# User Parameters
        //## Correspondences
        bool _symmetric = true;
        size_t _numNeighbours = 3;
        float _flagThreshold = 0.9f;
        //## Inliers
        float _kappaa = 3.0;
        bool _inlierUseOrientation = true;
        //## Transformation
        size_t _numIterations = 10;

        //# Internal Data structures

        //# Internal Parameters

        //# Internal functions
};

}//namespace registration

#endif // RIGIDREGISTRATION_HPP
