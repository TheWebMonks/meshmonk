#ifndef NONRIGIDREGISTRATION_HPP
#define NONRIGIDREGISTRATION_HPP

#include <Eigen/Dense>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <time.h>
#include <boost/variant.hpp>
#include "../global.hpp"
#include <CorrespondenceFilter.hpp>
#include <SymmetricCorrespondenceFilter.hpp>
#include <InlierDetector.hpp>
#include <ViscoElasticTransformer.hpp>

typedef Eigen::VectorXf VecDynFloat;
typedef Eigen::Matrix< float, Eigen::Dynamic, registration::NUM_FEATURES> FeatureMat; //matrix Mx6 of type float
typedef Eigen::Matrix< int, Eigen::Dynamic, 3> FacesMat;

namespace registration{

class NonrigidRegistration
{
    /*
    # GOAL
    This class performs icp-based nonrigid registration between two oriented pointclouds.

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
                       const FacesMat * const inFloatingFaces,
                       const VecDynFloat * const inFloatingFlags,
                       const VecDynFloat * const inTargetFlags);
        void set_parameters(bool symmetric,
                            size_t numNeighbours,
                            float kappaa,
                            size_t numIterations,
                            float sigmaSmoothing,
                            size_t numViscousIterationsStart,
                            size_t numViscousIterationsEnd,
                            size_t numElasticIterationsStart,
                            size_t numElasticIterationsEnd);

        void update();

    protected:

    private:
        //# Inputs/Outputs
        FeatureMat * _ioFloatingFeatures = NULL;
        const FeatureMat * _inTargetFeatures = NULL;
        const FacesMat * _inFloatingFaces;
        const VecDynFloat * _inFloatingFlags = NULL;
        const VecDynFloat * _inTargetFlags = NULL;

        //# User Parameters
        //## Correspondences
        bool _symmetric = true;
        size_t _numNeighbours = 3;
        //## Inliers
        float _kappaa = 3.0;
        //## Transformation
        size_t _numIterations = 10;
        float _sigmaSmoothing = 3.0;
        size_t _numViscousIterationsStart = 100;
        size_t _numViscousIterationsEnd = 1;
        size_t _numElasticIterationsStart = 100;
        size_t _numElasticIterationsEnd = 1;

        //# Internal Data structures

        //# Internal Parameters
        //## Transformation
        float _viscousAnnealingRate = exp(log(float(_numViscousIterationsEnd)/float(_numViscousIterationsStart))/_numIterations);
        float _elasticAnnealingRate = exp(log(float(_numElasticIterationsEnd)/float(_numElasticIterationsStart))/_numIterations);

        //# Internal functions
};

}//namespace registration

#endif // NONRIGIDREGISTRATION_HPP
