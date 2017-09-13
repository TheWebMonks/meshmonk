#ifndef PYRAMIDNONRIGIDREGISTRATION_HPP
#define PYRAMIDNONRIGIDREGISTRATION_HPP

#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <time.h>
#include <Eigen/Dense>
#include "../global.hpp"
#include "NonrigidRegistration.hpp"
#include "Downsampler.hpp"
#include "ScaleShifter.hpp"

typedef Eigen::VectorXf VecDynFloat;
typedef Eigen::VectorXi VecDynInt;
typedef Eigen::Matrix< float, Eigen::Dynamic, registration::NUM_FEATURES> FeatureMat; //matrix Mx6 of type float
typedef Eigen::Matrix< int, Eigen::Dynamic, 3> FacesMat;

namespace registration{

class PyramidNonrigidRegistration
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

        void set_input(FeatureMat &ioFloatingFeatures,
                       const FeatureMat &inTargetFeatures,
                       const FacesMat &inFloatingFaces,
                       const FacesMat &inTargetFaces,
                       const VecDynFloat &inFloatingFlags,
                       const VecDynFloat &inTargetFlags);

        void set_parameters(size_t numIterations = 60,
                            size_t numPyramidLayers = 3,
                            float downsampleFloatStart = 90.0f,
                            float downsampleTargetStart = 90.0f,
                            float downsampleFloatEnd = 0.0f,
                            float downsampleTargetEnd = 0.0f,
                            bool correspondencesSymmetric = true,
                            size_t correspondencesNumNeighbours = 5,
                            float _correspondencesFlagThreshold = 0.9f,
                            bool correspondencesEqualizePushPull = false,
                            float inlierKappa = 4.0f,
                            bool inlierUseOrientation = true,
                            float transformSigma = 3.0f,
                            size_t transformNumViscousIterationsStart = 200,
                            size_t transformNumViscousIterationsEnd = 1,
                            size_t transformNumElasticIterationsStart = 200,
                            size_t transformNumElasticIterationsEnd = 1);

        void update();

    protected:

    private:
        //# Inputs/Outputs
        FeatureMat * _ioFloatingFeatures = NULL;
        const FeatureMat * _inTargetFeatures = NULL;
        const FacesMat * _inFloatingFaces;
        const FacesMat * _inTargetFaces;
        const VecDynFloat * _inFloatingFlags = NULL;
        const VecDynFloat * _inTargetFlags = NULL;

        //# User Parameters
        //## Correspondences
        size_t _numIterations = 60;
        size_t _numPyramidLayers = 3;
        float _downsampleFloatStart = 90.0f; //percentage
        float _downsampleTargetStart = 90.0f; //percentage
        float _downsampleFloatEnd = 0.0f; //percentage
        float _downsampleTargetEnd = 0.0f; //percentage
        bool _correspondencesSymmetric = true;
        size_t _correspondencesNumNeighbours = 5;
        float _correspondencesFlagThreshold = 0.9f;
        bool _correspondencesEqualizePushPull = false;
        float _inlierKappa = 4.0f;
        bool _inlierUseOrientation = true;
        float _transformSigma = 3.0f;
        size_t _transformNumViscousIterationsStart = 200;
        size_t _transformNumViscousIterationsEnd = 1;
        size_t _transformNumElasticIterationsStart = 200;
        size_t _transformNumElasticIterationsEnd = 1;

        //# Internal Data structures

        //# Internal Parameters
        int _iterationsPerLayer = 0;
        //## Transformation
        float _viscousAnnealingRate = exp(log(float(_transformNumViscousIterationsEnd)/float(_transformNumViscousIterationsStart))/_numIterations);
        float _elasticAnnealingRate = exp(log(float(_transformNumElasticIterationsEnd)/float(_transformNumElasticIterationsStart))/_numIterations);
        std::vector<int> _viscousIterationsIntervals;
        std::vector<int> _elasticIterationsIntervals;
        //# Internal functions
};

}//namespace registration

#endif // PYRAMIDNONRIGIDREGISTRATION_HPP
