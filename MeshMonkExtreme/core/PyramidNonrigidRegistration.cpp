#include "PyramidPyramidNonrigidRegistration.hpp"

namespace registration {

void PyramidNonrigidRegistration::set_input(FeatureMat * const ioFloatingFeatures,
                            const FeatureMat * const inTargetFeatures,
                            const FacesMat * const inFloatingFaces,
                            const VecDynFloat * const inFloatingFlags,
                            const VecDynFloat * const inTargetFlags){
    _ioFloatingFeatures = ioFloatingFeatures;
    _inTargetFeatures = inTargetFeatures;
    _inFloatingFaces = inFloatingFaces;
    _inFloatingFlags = inFloatingFlags;
    _inTargetFlags = inTargetFlags;
}//end set_input()

void PyramidNonrigidRegistration::set_parameters(bool symmetric,
                            size_t numNeighbours,
                            float kappaa,
                            size_t numIterations,
                            float sigmaSmoothing,
                            size_t numViscousIterationsStart,
                            size_t numViscousIterationsEnd,
                            size_t numElasticIterationsStart,
                            size_t numElasticIterationsEnd){
    _symmetric = symmetric;
    _numNeighbours = numNeighbours;
    _kappaa = kappaa;
    _numIterations = numIterations;
    _sigmaSmoothing = sigmaSmoothing;
    _numViscousIterationsStart = numViscousIterationsStart;
    _numViscousIterationsEnd = numViscousIterationsEnd;
    _numElasticIterationsStart = numElasticIterationsStart;
    _numElasticIterationsEnd = numElasticIterationsEnd;
    _numViscousIterations = _numViscousIterationsStart;
    _numElasticIterations = _numElasticIterationsStart;

    _viscousAnnealingRate = exp(log(float(_numViscousIterationsEnd)/float(_numViscousIterationsStart))/(_numIterations-1));
    _elasticAnnealingRate = exp(log(float(_numElasticIterationsEnd)/float(_numElasticIterationsStart))/(_numIterations-1));
    //DEBUG
    std::cout << "viscous rate : " << _viscousAnnealingRate
    << " | start: " << _numViscousIterationsStart
    << " | end: " << _numViscousIterationsEnd
    << " | its: " << _numIterations << std::endl;
    std::cout << "elastic rate : " << _elasticAnnealingRate
    << " | start: " << _numElasticIterationsStart
    << " | end: " << _numElasticIterationsEnd
    << " | its: " << _numIterations << std::endl;
}//end set_parameters()


void PyramidNonrigidRegistration::update(){
    int numPyramidLayers = 3;
    int floatPyramidStart = 90;//percentage you want to start downsampling the floating mesh by
    int targetPyramidStart = 90;//percentage you want to start downsampling the target mesh by
    int floatPyramidEnd = 0;
    int targetPyramidEnd = 0;
    size_t numNonrigidIterations = 60;
    int numViscousIterationsStart = 500;
    int numElasticIterationsStart = 500;
    int numViscousIterationsStop = 1;
    int numElasticIterationsStop = 10;
    int iterationsPerLayer = int(std::round(float(numNonrigidIterations)/float(numPyramidLayers)));
    float viscousAnnealingRate = exp(log(float(numViscousIterationsStop)/float(numViscousIterationsStart))/numNonrigidIterations);
    float elasticAnnealingRate = exp(log(float(numElasticIterationsStop)/float(numElasticIterationsStart))/numNonrigidIterations);
    int viscousIterationsIntervals[numPyramidLayers + 1] = {};
    int elasticIterationsIntervals[numPyramidLayers + 1] = {};
    std::cout << "num viscous / elastic iterations : " << std::endl;
    for (int i = 0 ; i < numPyramidLayers ; i++) {
        viscousIterationsIntervals[i] = std::round(numViscousIterationsStart * pow(viscousAnnealingRate, i * iterationsPerLayer));
        elasticIterationsIntervals[i] = std::round(numElasticIterationsStart * pow(elasticAnnealingRate, i * iterationsPerLayer));
        std::cout << "num viscous iterations : " << viscousIterationsIntervals[i] << std::endl;
        std::cout << "num elastic iterations : " << elasticIterationsIntervals[i] << std::endl;
    }
    viscousIterationsIntervals[numPyramidLayers] = numViscousIterationsStop;
    elasticIterationsIntervals[numPyramidLayers] = numElasticIterationsStop;
    std::cout << "num viscous iterations : " << viscousIterationsIntervals[numPyramidLayers] << std::endl;
    std::cout << "num elastic iterations : " << elasticIterationsIntervals[numPyramidLayers] << std::endl;

    //# Initialize the floating features, their original indices and the faces.
    /*
    We need to do this before the pyramid iterations, because those have to be passed from the
    previous iteration to the next one. The 'ScaleShifter' class makes sure that the properties
    of the floating mesh of the previous pyramid scale are transferred to the current pyramid
    scale.
    */
    FeatureMat floatingFeatures;
    FacesMat floatingFaces;
    VecDynInt floatingOriginalIndices;
    for (int i = 0 ; i < numPyramidLayers ; i++){
        //# Copy the floating features and indices of the previous pyramid scale
        FeatureMat oldFloatingFeatures;
        VecDynInt oldFloatingOriginalIndices;
        if (i > 0) {
            oldFloatingFeatures = FeatureMat(floatingFeatures);
            oldFloatingOriginalIndices = VecDynInt(floatingOriginalIndices);
        }

        //# Downsample Floating Mesh
        float downsampleRatio = floatPyramidStart;
        if (numPyramidLayers > 1) {
            downsampleRatio = float(std::round(floatPyramidStart - i * std::round((floatPyramidStart-floatPyramidEnd)/(numPyramidLayers-1.0))));
        }
        downsampleRatio /= 100.0f;
        std::cout<< " DOWNSAMPLE RATIO: " << downsampleRatio << std::endl;
        registration::Downsampler downsampler;
        VecDynFloat floatingFlags;
        downsampler.set_input(&originalFloatingFeatures, &originalFloatingFaces, &originalFloatingFlags);
        downsampler.set_output(floatingFeatures, floatingFaces, floatingFlags, floatingOriginalIndices);
        downsampler.set_parameters(downsampleRatio);
        downsampler.update();

        //# Downsample Target Mesh
        downsampleRatio = targetPyramidStart;
        if (numPyramidLayers > 1){
            downsampleRatio = float(std::round(targetPyramidStart - i * std::round((targetPyramidStart-targetPyramidEnd)/(numPyramidLayers-1.0))));
        }
        downsampleRatio /= 100.0f;
        FeatureMat targetFeatures;
        FacesMat targetFaces;
        VecDynFloat targetFlags;
        downsampler.set_input(&originalTargetFeatures, &originalTargetFaces, &originalTargetFlags);
        downsampler.set_output(targetFeatures, targetFaces, targetFlags);
        downsampler.set_parameters(downsampleRatio);
        downsampler.update();

        //# Transfer floating mesh properties from previous pyramid scale to the current one.
        if (i > 0) {
            //## Scale up
            registration::ScaleShifter scaleShifter;
            scaleShifter.set_input(oldFloatingFeatures, oldFloatingOriginalIndices, floatingOriginalIndices);
            scaleShifter.set_output(floatingFeatures);
            scaleShifter.update();
        }

        //# Registration
        //DEBUG
        float sigmaSmoothing = 10.0f;
        size_t numTargetVertices = targetFeatures.rows();
        registration::NonrigidRegistration nonrigidRegistration;
        nonrigidRegistration.set_input(&floatingFeatures, &targetFeatures, &floatingFaces, &floatingFlags, &targetFlags);
        nonrigidRegistration.set_parameters(true, 5, 3.0,
                                            iterationsPerLayer, sigmaSmoothing,
                                            viscousIterationsIntervals[i], viscousIterationsIntervals[i+1],
                                            elasticIterationsIntervals[i], elasticIterationsIntervals[i+1]);
        nonrigidRegistration.update();
    }

}//end update()

}//namespace registration
