#include "PyramidNonrigidRegistration.hpp"

namespace registration {

void PyramidNonrigidRegistration::set_input(FeatureMat &ioFloatingFeatures,
                                       const FeatureMat &inTargetFeatures,
                                       const FacesMat &inFloatingFaces,
                                       const FacesMat &inTargetFaces,
                                       const VecDynFloat &inFloatingFlags,
                                       const VecDynFloat &inTargetFlags){
    _ioFloatingFeatures = &ioFloatingFeatures;
    _inTargetFeatures = &inTargetFeatures;
    _inFloatingFaces = &inFloatingFaces;
    _inTargetFaces = &inTargetFaces;
    _inFloatingFlags = &inFloatingFlags;
    _inTargetFlags = &inTargetFlags;
}//end set_input()

void PyramidNonrigidRegistration::set_parameters(size_t numIterations /*= 60*/,
                                                size_t numPyramidLayers /* = 3*/,
                                                float downsampleFloatStart /* = 90.0f*/,
                                                float downsampleTargetStart /* = 90.0f*/,
                                                float downsampleFloatEnd /* = 0.0f*/,
                                                float downsampleTargetEnd /* = 0.0f*/,
                                                bool correspondencesSymmetric /* = true*/,
                                                size_t correspondencesNumNeighbours /* = 5*/,
                                                float inlierKappa /* = 4.0f*/,
                                                bool inlierUseOrientation /* = true*/,
                                                float transformSigma /* = 3.0f*/,
                                                size_t transformNumViscousIterationsStart /* = 200*/,
                                                size_t transformNumViscousIterationsEnd /* = 1*/,
                                                size_t transformNumElasticIterationsStart /* = 200*/,
                                                size_t transformNumElasticIterationsEnd /* = 1*/){

        //# User Parameters
        _numIterations = numIterations;
        if (_numIterations <= 0) {
            _numIterations = 1;
            std::cout << "Number of iterations has to be a positive integer larger than 0!" <<std::endl;
        }
        _numPyramidLayers = numPyramidLayers;
        if (_numPyramidLayers <= 0) {
            _numPyramidLayers = 1;
            std::cout << "Number of pyramid layers has to be a positive integer larger than 0!" <<std::endl;
        }
        _downsampleFloatStart = downsampleFloatStart; //percentage
        if ((_downsampleFloatStart < 0.0f) || (_downsampleFloatStart >= 100.0f)) {
            _downsampleFloatStart = 90.0f;
            std::cout << "Downsample percentages have to be larger or equal to 0.0 and smaller than 100.0" <<std::endl;
        }
        _downsampleTargetStart = downsampleTargetStart; //percentage
        if ((_downsampleTargetStart < 0.0f) || (_downsampleFloatStart >= 100.0f)) {
            _downsampleTargetStart = 90.0f;
            std::cout << "Downsample percentages have to be larger or equal to 0.0 and smaller than 100.0" <<std::endl;
        }
        _downsampleFloatEnd = downsampleFloatEnd; //percentage
        _downsampleTargetEnd = downsampleTargetEnd; //percentage
        _correspondencesSymmetric = correspondencesSymmetric;
        _correspondencesNumNeighbours = correspondencesNumNeighbours;
        _inlierKappa = inlierKappa;
        _inlierUseOrientation = inlierUseOrientation;
        _transformSigma = transformSigma;
        _transformNumViscousIterationsStart = transformNumViscousIterationsStart;
        _transformNumViscousIterationsEnd = transformNumViscousIterationsEnd;
        _transformNumElasticIterationsStart = transformNumElasticIterationsStart;
        _transformNumElasticIterationsEnd = transformNumElasticIterationsEnd;

        //# Internal Parameters
        //## Determine annealing rates
        _iterationsPerLayer = int(std::round(float(_numIterations)/float(_numPyramidLayers)));
        if (_transformNumViscousIterationsStart > 0){
            _viscousAnnealingRate = exp(log((float(_transformNumViscousIterationsEnd)+0.5f)/float(_transformNumViscousIterationsStart))/_numIterations);//the +0.5f is for numerical stability
        }
        else {_viscousAnnealingRate = 1.0f;}
        if (_transformNumElasticIterationsStart > 0){
        _elasticAnnealingRate = exp(log((float(_transformNumElasticIterationsEnd)+0.5f)/float(_transformNumElasticIterationsStart))/_numIterations);//the +0.5f is for numerical stability
        }
        else {_elasticAnnealingRate = 1.0f;}

        //## Determine smoothing iterations for each pyramid layer
        _viscousIterationsIntervals.resize(_numPyramidLayers + 1);
        _elasticIterationsIntervals.resize(_numPyramidLayers + 1);
        std::cout << "viscous / elastic annealing rate : " << _viscousAnnealingRate << " / " << _elasticAnnealingRate << std::endl;
        std::cout << "num viscous / elastic iterations : " << std::endl;
        for (size_t i = 0 ; i < _numPyramidLayers ; i++) {
            _viscousIterationsIntervals[i] = std::round(_transformNumViscousIterationsStart * pow(_viscousAnnealingRate, i * _iterationsPerLayer));
            _elasticIterationsIntervals[i] = std::round(_transformNumElasticIterationsStart * pow(_elasticAnnealingRate, i * _iterationsPerLayer));
            std::cout << "num viscous iterations : " << _viscousIterationsIntervals[i] << std::endl;
            std::cout << "num elastic iterations : " << _elasticIterationsIntervals[i] << std::endl;
        }
        _viscousIterationsIntervals[_numPyramidLayers] = _transformNumViscousIterationsEnd;
        _elasticIterationsIntervals[_numPyramidLayers] = _transformNumElasticIterationsEnd;
        std::cout << "num viscous iterations : " << _viscousIterationsIntervals[_numPyramidLayers] << std::endl;
        std::cout << "num elastic iterations : " << _elasticIterationsIntervals[_numPyramidLayers] << std::endl;
}//end set_parameters()


void PyramidNonrigidRegistration::update(){

    //# Initialize the floating features, their original indices and the faces.
    /*
    We need to do this before the pyramid iterations, because those have to be passed from the
    previous iteration to the next one. The 'ScaleShifter' class makes sure that the properties
    of the floating mesh of the previous pyramid scale are transferred to the current pyramid
    scale.
    */
    size_t numFloatingFeatures = _ioFloatingFeatures->rows();
    FeatureMat floatingFeatures;
    FacesMat floatingFaces;
    VecDynInt floatingOriginalIndices;
    FeatureMat oldFloatingFeatures;
    VecDynInt oldFloatingOriginalIndices;

    //# Start Pyramid Nonrigid Registration
    for (size_t i = 0 ; i < _numPyramidLayers ; i++){

        //# Downsample Floating Mesh
        //## Determine the downsample ratio for the current pyramid layer
        float downsampleRatio = _downsampleFloatStart;
        if (_numPyramidLayers > 1) {
            downsampleRatio = float(std::round(_downsampleFloatStart - i * std::round((_downsampleFloatStart-_downsampleFloatEnd)/(_numPyramidLayers-1.0))));
        }
        downsampleRatio /= 100.0f;
        std::cout<< " DOWNSAMPLE RATIO       : " << downsampleRatio << std::endl;
        //## Set up Downsampler
        Downsampler downsampler;
        VecDynFloat floatingFlags;
        downsampler.set_input(_ioFloatingFeatures, _inFloatingFaces, _inFloatingFlags);
        downsampler.set_output(floatingFeatures, floatingFaces, floatingFlags, floatingOriginalIndices);
        downsampler.set_parameters(downsampleRatio);
        downsampler.update();

        //# Downsample Target Mesh
        //## Determine the downsample ratio for the current pyramid layer
        downsampleRatio = _downsampleTargetStart;
        if (_numPyramidLayers > 1){
            downsampleRatio = float(std::round(_downsampleTargetStart - i * std::round((_downsampleTargetStart-_downsampleTargetEnd)/(_numPyramidLayers-1.0))));
        }
        downsampleRatio /= 100.0f;
        //## Set up Downsampler
        FeatureMat targetFeatures;
        FacesMat targetFaces;
        VecDynFloat targetFlags;
        downsampler.set_input(_inTargetFeatures, _inTargetFaces, _inTargetFlags);
        downsampler.set_output(targetFeatures, targetFaces, targetFlags);
        downsampler.set_parameters(downsampleRatio);
        downsampler.update();

        //# Transfer floating mesh properties from previous pyramid scale to the current one.
        if (i > 0) {
            //## Scale up
            ScaleShifter scaleShifter;
            scaleShifter.set_input(oldFloatingFeatures, oldFloatingOriginalIndices, floatingOriginalIndices);
            scaleShifter.set_output(floatingFeatures);
            scaleShifter.update();
        }

        //# Registration
        NonrigidRegistration nonrigidRegistration;
        nonrigidRegistration.set_input(&floatingFeatures, &targetFeatures, &floatingFaces, &floatingFlags, &targetFlags);
        nonrigidRegistration.set_parameters(_correspondencesSymmetric, _correspondencesNumNeighbours,
                                            _inlierKappa, _inlierUseOrientation,
                                            _iterationsPerLayer, _transformSigma,
                                            _viscousIterationsIntervals[i], _viscousIterationsIntervals[i+1],
                                            _elasticIterationsIntervals[i], _elasticIterationsIntervals[i+1]);
        nonrigidRegistration.update();

        //# Copy the result in temporary variables for use in the next pyramid scale
        oldFloatingFeatures = FeatureMat(floatingFeatures);
        oldFloatingOriginalIndices = VecDynInt(floatingOriginalIndices);
    }// Pyramid iteratations

    std::cout << "A - _ioFloatingFeatures->topRows(10) : \n" << _ioFloatingFeatures->topRows(10) << std::endl;
    //# Copy result to output
    VecDynInt originalIndices = VecDynInt::Zero(numFloatingFeatures);
    for (size_t j = 0 ; j < numFloatingFeatures ; j++){ originalIndices(j) = j; }
    ScaleShifter scaleShifter;
    scaleShifter.set_input(floatingFeatures, floatingOriginalIndices, originalIndices);
    scaleShifter.set_output(*_ioFloatingFeatures);
    scaleShifter.update();

    std::cout << "B - _ioFloatingFeatures->topRows(10) : \n" << _ioFloatingFeatures->topRows(10) << std::endl;

}//end update()

}//namespace registration
