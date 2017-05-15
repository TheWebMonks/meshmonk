#include "NonrigidRegistration.hpp"

namespace registration {

void NonrigidRegistration::set_input(FeatureMat * const ioFloatingFeatures,
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

void NonrigidRegistration::set_parameters(bool symmetric,
                            size_t numNeighbours,
                            float flagThreshold,
                            float kappaa,
                            bool inlierUseOrientation,
                            size_t numIterations,
                            float sigmaSmoothing,
                            size_t numViscousIterationsStart,
                            size_t numViscousIterationsEnd,
                            size_t numElasticIterationsStart,
                            size_t numElasticIterationsEnd){
    _symmetric = symmetric;
    _numNeighbours = numNeighbours;
    _flagThreshold = flagThreshold;
    _kappaa = kappaa;
    _inlierUseOrientation = inlierUseOrientation;
    _numIterations = numIterations;
    _sigmaSmoothing = sigmaSmoothing;
    _numViscousIterationsStart = numViscousIterationsStart;
    _numViscousIterationsEnd = numViscousIterationsEnd;
    _numElasticIterationsStart = numElasticIterationsStart;
    _numElasticIterationsEnd = numElasticIterationsEnd;
    _numViscousIterations = _numViscousIterationsStart;
    _numElasticIterations = _numElasticIterationsStart;

    //DEBUG
        std::cout << "NonrigidRegistration::set_parameters parameters: \n"
        << " numViscousIterationsStart: " << numViscousIterationsStart << "\n"
        << " numViscousIterationsEnd: " << numViscousIterationsEnd << "\n"
        << " numElasticIterationsStart: " << numElasticIterationsStart << "\n"
        << " numElasticIterationsEnd: " << numElasticIterationsEnd << std::endl;
        //END DEBUG

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


void NonrigidRegistration::update(){

    //# Initializes
    size_t numFloatingVertices = _ioFloatingFeatures->rows();
    FeatureMat correspondingFeatures = FeatureMat::Zero(numFloatingVertices, registration::NUM_FEATURES);
    VecDynFloat correspondingFlags = VecDynFloat::Zero(numFloatingVertices);

    //# Set up the filters
    //## Correspondence Filter (note: this part should be improved by using inheritance in the correspondence filter classes)

    BaseCorrespondenceFilter* correspondenceFilter = NULL;
    //std::unique_ptr<BaseCorrespondenceFilter> correspondenceFilter = new SymmetricCorrespondenceFilter;

    if (_symmetric) {
        correspondenceFilter = new SymmetricCorrespondenceFilter();
    }
    else {
        correspondenceFilter = new CorrespondenceFilter();
    }
    correspondenceFilter->set_floating_input(_ioFloatingFeatures, _inFloatingFlags);
    correspondenceFilter->set_target_input(_inTargetFeatures, _inTargetFlags);
    correspondenceFilter->set_output(&correspondingFeatures, &correspondingFlags);
    correspondenceFilter->set_parameters(_numNeighbours, _flagThreshold);

    //## Inlier Filter
    VecDynFloat floatingWeights = VecDynFloat::Ones(numFloatingVertices);
    InlierDetector inlierDetector;
    inlierDetector.set_input(_ioFloatingFeatures, &correspondingFeatures,
                                &correspondingFlags);
    inlierDetector.set_output(&floatingWeights);
    inlierDetector.set_parameters(_kappaa, _inlierUseOrientation);
    //## Transformation Filter
    _numViscousIterations = _numViscousIterationsStart;
    _numElasticIterations = _numElasticIterationsStart;
    ViscoElasticTransformer transformer;
    transformer.set_input(&correspondingFeatures, &floatingWeights, _inFloatingFlags, _inFloatingFaces);
    transformer.set_output(_ioFloatingFeatures);

    //# Perform ICP
    time_t timeStart, timePreIteration, timePostIteration, timeEnd;
    timeStart = time(0);
    std::cout << "Starting Nonrigid Registration process..." << std::endl;
    for (size_t iteration = 0 ; iteration < _numIterations ; iteration++) {
        timePreIteration = time(0);

        //# Anneal parameters
        _numViscousIterations = int(std::round(_numViscousIterationsStart * std::pow(_viscousAnnealingRate, iteration)));
        if (_numViscousIterations < _numViscousIterationsEnd) { _numViscousIterations = _numViscousIterationsEnd;}
        _numElasticIterations = int(std::round(_numElasticIterationsStart * std::pow(_elasticAnnealingRate, iteration)));
        if (_numElasticIterations < _numElasticIterationsEnd) { _numElasticIterations = _numElasticIterationsEnd;}


        //# Correspondences
        correspondenceFilter->set_floating_input(_ioFloatingFeatures, _inFloatingFlags);
        correspondenceFilter->set_target_input(_inTargetFeatures, _inTargetFlags);
        correspondenceFilter->update();

        //# Inlier Detection
        inlierDetector.update();

        //# Transformation
        transformer.set_parameters(10, _sigmaSmoothing, _numViscousIterations,_numElasticIterations);
        transformer.update();

        //# Print info
        timePostIteration = time(0);
        std::cout << "Iteration " << iteration+1 << "/" << _numIterations << " took "<< difftime(timePostIteration, timePreIteration) <<" second(s)."<< std::endl;
    }
    timeEnd = time(0);
    std::cout << "Nonrigid Registration Completed in " << difftime(timeEnd, timeStart) <<" second(s)."<< std::endl;

    delete correspondenceFilter;

}//end update()

}//namespace registration
