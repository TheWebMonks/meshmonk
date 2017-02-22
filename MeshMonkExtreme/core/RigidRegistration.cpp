#include "RigidRegistration.hpp"

namespace registration {

void RigidRegistration::update(){

    //# Initializes
    size_t numFloatingVertices = _ioFloatingFeatures->rows();
    size_t numTargetVertices = _inTargetFeatures->rows();
    FeatureMat correspondingFeatures = FeatureMat::Zero(numFloatingVertices, registration::NUM_FEATURES);
    VecDynFloat correspondingFlags = VecDynFloat::Zero(numFloatingVertices);

    //# Set up the filter
    //## Correspondence Filter (note: this part should be improved by using inheritance in the correspondence filter classes)
    BaseCorrespondenceFilter correspondenceFilter;
        if (_symmetric) {
        correspondenceFilter = SymmetricCorrespondenceFilter();
    }
    else {
        correspondenceFilter = CorrespondenceFilter();
    }
    correspondenceFilter.set_floating_input(_ioFloatingFeatures, _inFloatingFlags);
    correspondenceFilter.set_target_input(_inTargetFeatures, _inTargetFlags);
    correspondenceFilter.set_output(&correspondingFeatures, &correspondingFlags);
    correspondenceFilter.set_parameters(_numNeighbours);

    //## Inlier Filter
    VecDynFloat floatingWeights = VecDynFloat::Ones(numFloatingVertices);
    InlierDetector inlierDetector;
    inlierDetector.set_input(_ioFloatingFeatures, &correspondingFeatures,
                                &correspondingFlags);
    inlierDetector.set_output(&floatingWeights);
    inlierDetector.set_parameters(3.0f);
    //## Transformation Filter
    RigidTransformer rigidTransformer;
    rigidTransformer.set_input(&correspondingFeatures, &floatingWeights);
    rigidTransformer.set_output(_ioFloatingFeatures);
    rigidTransformer.set_parameters(false);

//## Initialize
//        numFloatingVertices = self.ioFloatingFeatures.shape[0]
//        numTargetVertices = self.inTargetFeatures.shape[0]
//        floatingFlags = numpy.ones((numFloatingVertices), dtype = float)
//        targetFlags = numpy.ones((numTargetVertices), dtype = float)
//        correspondingFeatures = numpy.zeros((self.ioFloatingFeatures.shape), dtype = float)
//        correspondingFlags = numpy.ones((floatingFlags.shape), dtype = float)
//        correspondenceFilter = []
//        if (self.symmetric):
//            correspondenceFilter = SymCorrespondenceFilter(self.ioFloatingFeatures,
//                                                           floatingFlags,
//                                                           self.inTargetFeatures,
//                                                           targetFlags,
//                                                           correspondingFeatures,
//                                                           correspondingFlags,
//                                                           self.wknnNumNeighbours,
//                                                           self.ratioPushToPull)
//        else:
//            correspondenceFilter = CorrespondenceFilter(self.ioFloatingFeatures,
//                                                        floatingFlags,
//                                                        self.inTargetFeatures,
//                                                        targetFlags,
//                                                        correspondingFeatures,
//                                                        correspondingFlags,
//                                                        self.wknnNumNeighbours)
//
//        ## Set up inlier filter
//        floatingWeights = numpy.ones((numFloatingVertices), dtype = float)
//        inlierFilter = InlierFilter(self.ioFloatingFeatures, correspondingFeatures,
//                                    correspondingFlags, floatingWeights,
//                                    self.kappaa)
//
//
//        ## Set up transformation filter
//        transformationFilter = RigidTransformationFilter(self.ioFloatingFeatures,
//                                                         correspondingFeatures,
//                                                         floatingWeights)
//
//
//        print ("Starting Rigid Registration...")
//        timePre = time.time()
//        for iteration in range(self.numIterations):
//            timeStart = time.time()
//            ## 1) Update Nearest neighbours.
//            correspondenceFilter.set_floating_features(self.ioFloatingFeatures, floatingFlags)
//            correspondenceFilter.update()
//            ## 2) Update inlier weights.
//            inlierFilter.update()
//            ## 3) Update transformation.
//            transformationFilter.update()
//
//            ## Print progress
//            timeEnd = time.time()
//            print("Iteration " + str(iteration) + "/" + str(self.numIterations) + " took " + str(timeEnd-timeStart) + 's')
//            iteration = iteration + 1
//
//        timePost = time.time()
//        print ("Completed Rigid Registration in " + str(timePost-timePre) + 's')
}

}//namespace registration
