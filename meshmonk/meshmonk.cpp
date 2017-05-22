#include "meshmonk.hpp"

namespace meshmonk{

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus


    //######################################################################################
    //################################  MEX WRAPPING  ######################################
    //######################################################################################
    /*
    We're wrapping some functionality in the library so it can be easily mexed in Matlab
    */

    void pyramid_registration_mex(float floatingFeaturesArray[], const float targetFeaturesArray[],
                                const size_t numFloatingElements, const size_t numTargetElements,
                                const int floatingFacesArray[], const int targetFacesArray[],
                                const size_t numFloatingFaces, const size_t numTargetFaces,
                                const float floatingFlagsArray[], const float targetFlagsArray[],
                                const size_t numIterations/*= 60*/, const size_t numPyramidLayers/*= 3*/,
                                const float downsampleFloatStart/*= 90*/, const float downsampleTargetStart/*= 90*/,
                                const float downsampleFloatEnd/*= 0*/, const float downsampleTargetEnd/*= 0*/,
                                const bool correspondencesSymmetric/*= true*/, const size_t correspondencesNumNeighbours/*= 5*/,
                                const float correspondencesFlagThreshold/* = 0.9f*/,
                                const float inlierKappa/*= 4.0f*/, const bool inlierUseOrientation/*=true*/,
                                const float transformSigma/*= 3.0f*/,
                                const size_t transformNumViscousIterationsStart/*= 50*/, const size_t transformNumViscousIterationsEnd/*= 1*/,
                                const size_t transformNumElasticIterationsStart/*= 50*/, const size_t transformNumElasticIterationsEnd/*= 1*/){
        //# Convert arrays to Eigen matrices (see http://dovgalecs.com/blog/eigen-how-to-get-in-and-out-data-from-eigen-matrix/)
        FeatureMat floatingFeatures = Eigen::Map<FeatureMat>(floatingFeaturesArray, numFloatingElements, registration::NUM_FEATURES);
        const FeatureMat targetFeatures = Eigen::Map<const FeatureMat>(targetFeaturesArray, numTargetElements, registration::NUM_FEATURES);
        const FacesMat floatingFaces = Eigen::Map<const FacesMat>(floatingFacesArray, numFloatingFaces, 3);
        const FacesMat targetFaces = Eigen::Map<const FacesMat>(targetFacesArray, numTargetFaces, 3);
        const VecDynFloat floatingFlags = Eigen::Map<const VecDynFloat>(floatingFlagsArray, numFloatingElements);
        const VecDynFloat targetFlags = Eigen::Map<const VecDynFloat>(targetFlagsArray, numTargetElements);

        pyramid_registration(floatingFeatures, targetFeatures,
                                floatingFaces, targetFaces,
                                floatingFlags, targetFlags,
                                numIterations, numPyramidLayers,
                                downsampleFloatStart, downsampleTargetStart,
                                downsampleFloatEnd, downsampleTargetEnd,
                                correspondencesSymmetric, correspondencesNumNeighbours,
                                correspondencesFlagThreshold,
                                inlierKappa, inlierUseOrientation,
                                transformSigma,
                                transformNumViscousIterationsStart, transformNumViscousIterationsEnd,
                                transformNumElasticIterationsStart, transformNumElasticIterationsEnd);

        //# Convert back to raw data
        Eigen::Map<FeatureMat>(floatingFeaturesArray, floatingFeatures.rows(), floatingFeatures.cols()) = floatingFeatures;

    }


    void nonrigid_registration_mex(float floatingFeaturesArray[], const float targetFeaturesArray[],
                                const size_t numFloatingElements, const size_t numTargetElements,
                                const int floatingFacesArray[], const int targetFacesArray[],
                                const size_t numFloatingFaces, const size_t numTargetFaces,
                                const float floatingFlagsArray[], const float targetFlagsArray[],
                                const size_t numIterations/*= 60*/,
                                const bool correspondencesSymmetric/*= true*/, const size_t correspondencesNumNeighbours/*= 5*/,
                                const float correspondencesFlagThreshold/* = 0.9f*/,
                                const float inlierKappa/*= 4.0f*/, const bool inlierUseOrientation/*=true*/,
                                const float transformSigma/*= 3.0f*/,
                                const size_t transformNumViscousIterationsStart/*= 50*/, const size_t transformNumViscousIterationsEnd/*= 1*/,
                                const size_t transformNumElasticIterationsStart/*= 50*/, const size_t transformNumElasticIterationsEnd/*= 1*/){
        //# Convert arrays to Eigen matrices (see http://dovgalecs.com/blog/eigen-how-to-get-in-and-out-data-from-eigen-matrix/)
        FeatureMat floatingFeatures = Eigen::Map<FeatureMat>(floatingFeaturesArray, numFloatingElements, registration::NUM_FEATURES);
        const FeatureMat targetFeatures = Eigen::Map<const FeatureMat>(targetFeaturesArray, numTargetElements, registration::NUM_FEATURES);
        const FacesMat floatingFaces = Eigen::Map<const FacesMat>(floatingFacesArray, numFloatingFaces, 3);
        const FacesMat targetFaces = Eigen::Map<const FacesMat>(targetFacesArray, numTargetFaces, 3);
        const VecDynFloat floatingFlags = Eigen::Map<const VecDynFloat>(floatingFlagsArray, numFloatingElements);
        const VecDynFloat targetFlags = Eigen::Map<const VecDynFloat>(targetFlagsArray, numTargetElements);

        //# Run nonrigid registration
        nonrigid_registration(floatingFeatures, targetFeatures,
                                floatingFaces, targetFaces,
                                floatingFlags, targetFlags,
                                numIterations,
                                correspondencesSymmetric, correspondencesNumNeighbours,
                                correspondencesFlagThreshold,
                                inlierKappa, inlierUseOrientation,
                                transformSigma,
                                transformNumViscousIterationsStart, transformNumViscousIterationsEnd,
                                transformNumElasticIterationsStart, transformNumElasticIterationsEnd);

        //# Convert back to raw data
        Eigen::Map<FeatureMat>(floatingFeaturesArray, floatingFeatures.rows(), floatingFeatures.cols()) = floatingFeatures;
    }


    void rigid_registration_mex(float floatingFeaturesArray[], const float targetFeaturesArray[],
                                const size_t numFloatingElements, const size_t numTargetElements,
                                const int floatingFacesArray[], const int targetFacesArray[],
                                const size_t numFloatingFaces, const size_t numTargetFaces,
                                const float floatingFlagsArray[], const float targetFlagsArray[],
                                const size_t numIterations/*= 60*/,
                                const bool correspondencesSymmetric/*= true*/, const size_t correspondencesNumNeighbours/*= 5*/,
                                const float correspondencesFlagThreshold/* = 0.9f*/,
                                const float inlierKappa/*= 4.0f*/, const bool inlierUseOrientation/*=true*/){
        //# Convert arrays to Eigen matrices (see http://dovgalecs.com/blog/eigen-how-to-get-in-and-out-data-from-eigen-matrix/)
        FeatureMat floatingFeatures = Eigen::Map<FeatureMat>(floatingFeaturesArray, numFloatingElements, registration::NUM_FEATURES);
        const FeatureMat targetFeatures = Eigen::Map<const FeatureMat>(targetFeaturesArray, numTargetElements, registration::NUM_FEATURES);
        const FacesMat floatingFaces = Eigen::Map<const FacesMat>(floatingFacesArray, numFloatingFaces, 3);
        const FacesMat targetFaces = Eigen::Map<const FacesMat>(targetFacesArray, numTargetFaces, 3);
        const VecDynFloat floatingFlags = Eigen::Map<const VecDynFloat>(floatingFlagsArray, numFloatingElements);
        const VecDynFloat targetFlags = Eigen::Map<const VecDynFloat>(targetFlagsArray, numTargetElements);

        //# Run rigid registration
        rigid_registration(floatingFeatures, targetFeatures,
                            floatingFaces, targetFaces,
                            floatingFlags, targetFlags,
                            numIterations,
                            correspondencesSymmetric, correspondencesNumNeighbours,
                            correspondencesFlagThreshold,
                            inlierKappa);

        //# Convert back to raw data
        Eigen::Map<FeatureMat>(floatingFeaturesArray, floatingFeatures.rows(), floatingFeatures.cols()) = floatingFeatures;
    }


    void compute_correspondences_mex(const float floatingFeaturesArray[], const float targetFeaturesArray[],
                                    const size_t numFloatingElements, const size_t numTargetElements,
                                    const float floatingFlagsArray[], const float targetFlagsArray[],
                                    float correspondingFeaturesArray[], float correspondingFlagsArray[],
                                    const bool correspondencesSymmetric/*= true*/, const size_t correspondencesNumNeighbours/*= 5*/,
                                    const float correspondencesFlagThreshold /*= 0.9f*/){
        //# Convert arrays to Eigen matrices (see http://dovgalecs.com/blog/eigen-how-to-get-in-and-out-data-from-eigen-matrix/)
        const FeatureMat floatingFeatures = Eigen::Map<const FeatureMat>(floatingFeaturesArray, numFloatingElements, registration::NUM_FEATURES);
        const FeatureMat targetFeatures = Eigen::Map<const FeatureMat>(targetFeaturesArray, numTargetElements, registration::NUM_FEATURES);
        const VecDynFloat floatingFlags = Eigen::Map<const VecDynFloat>(floatingFlagsArray, numFloatingElements);
        const VecDynFloat targetFlags = Eigen::Map<const VecDynFloat>(targetFlagsArray, numTargetElements);
        FeatureMat correspondingFeatures = Eigen::Map<FeatureMat>(correspondingFeaturesArray, numFloatingElements, registration::NUM_FEATURES);
        VecDynFloat correspondingFlags = Eigen::Map<VecDynFloat>(correspondingFlagsArray, numFloatingElements);

        //# Compute Correspondences
        compute_correspondences(floatingFeatures, targetFeatures,
                                floatingFlags, targetFlags,
                                correspondingFeatures, correspondingFlags,
                                correspondencesSymmetric, correspondencesNumNeighbours,
                                correspondencesFlagThreshold);

        //# Convert back to raw data
        Eigen::Map<FeatureMat>(correspondingFeaturesArray, correspondingFeatures.rows(), correspondingFeatures.cols()) = correspondingFeatures;
        Eigen::Map<VecDynFloat>(correspondingFlagsArray, correspondingFlags.rows(), correspondingFlags.cols()) = correspondingFlags;
    }


    void compute_inlier_weights_mex(const float floatingFeaturesArray[], const float correspondingFeaturesArray[],
                                    const size_t numFloatingElements,
                                    const float correspondingFlagsArray[], float inlierWeightsArray[],
                                    const float inlierKappa/*= 4.0f*/, const bool useOrientation/*= true*/){
        //# Convert arrays to Eigen matrices (see http://dovgalecs.com/blog/eigen-how-to-get-in-and-out-data-from-eigen-matrix/)
        const FeatureMat floatingFeatures = Eigen::Map<const FeatureMat>(floatingFeaturesArray, numFloatingElements, registration::NUM_FEATURES);
        const FeatureMat correspondingFeatures = Eigen::Map<const FeatureMat>(correspondingFeaturesArray, numFloatingElements, registration::NUM_FEATURES);
        const VecDynFloat correspondingFlags = Eigen::Map<const VecDynFloat>(correspondingFlagsArray, numFloatingElements);
        VecDynFloat inlierWeights = Eigen::Map<VecDynFloat>(inlierWeightsArray, numFloatingElements);

        //# Computer Inlier Weights
        compute_inlier_weights(floatingFeatures, correspondingFeatures,
                                correspondingFlags, inlierWeights,
                                inlierKappa, useOrientation);

        //# Convert back to raw data
        Eigen::Map<VecDynFloat>(inlierWeightsArray, inlierWeights.rows(), inlierWeights.cols()) = inlierWeights;
    }


    void compute_rigid_transformation_mex(float floatingFeaturesArray[], const size_t numFloatingElements,
                                        const float correspondingFeaturesArray[], const float inlierWeightsArray[],
                                        const bool allowScaling /*= fakse*/){
        //# Convert arrays to Eigen matrices (see http://dovgalecs.com/blog/eigen-how-to-get-in-and-out-data-from-eigen-matrix/)
        FeatureMat floatingFeatures = Eigen::Map<FeatureMat>(floatingFeaturesArray, numFloatingElements, registration::NUM_FEATURES);
        const FeatureMat correspondingFeatures = Eigen::Map<const FeatureMat>(correspondingFeaturesArray, numFloatingElements, registration::NUM_FEATURES);
        const VecDynFloat inlierWeights = Eigen::Map<const VecDynFloat>(inlierWeightsArray, numFloatingElements);

        //# Run nonrigid registration
        compute_rigid_transformation(floatingFeatures, correspondingFeatures,
                                    inlierWeights, allowScaling);

        //# Convert back to raw data
        Eigen::Map<FeatureMat>(floatingFeaturesArray, floatingFeatures.rows(), floatingFeatures.cols()) = floatingFeatures;
    }


    void compute_nonrigid_transformation_mex(float floatingFeaturesArray[], const float correspondingFeaturesArray[],
                                            const size_t numFloatingElements,
                                            const int floatingFacesArray[], const size_t numFloatingFaces,
                                            const float floatingFlagsArray[], const float inlierWeightsArray[],
                                            const size_t transformNumNeighbours/*= 10*/, const float transformSigma/*= 3.0f*/,
                                            const size_t transformNumViscousIterations/*= 50*/, const size_t transformNumElasticIterations/*= 50*/){
        //# Convert arrays to Eigen matrices (see http://dovgalecs.com/blog/eigen-how-to-get-in-and-out-data-from-eigen-matrix/)
        FeatureMat floatingFeatures = Eigen::Map<FeatureMat>(floatingFeaturesArray, numFloatingElements, registration::NUM_FEATURES);
        const FeatureMat correspondingFeatures = Eigen::Map<const FeatureMat>(correspondingFeaturesArray, numFloatingElements, registration::NUM_FEATURES);
        const FacesMat floatingFaces = Eigen::Map<const FacesMat>(floatingFacesArray, numFloatingFaces, 3);
        const VecDynFloat floatingFlags = Eigen::Map<const VecDynFloat>(floatingFlagsArray, numFloatingElements);
        const VecDynFloat inlierWeights = Eigen::Map<const VecDynFloat>(inlierWeightsArray, numFloatingElements);

        //# Run nonrigid registration
        compute_nonrigid_transformation(floatingFeatures, correspondingFeatures,
                                        floatingFaces, floatingFlags,
                                        inlierWeights,
                                        transformNumNeighbours, transformSigma,
                                        transformNumViscousIterations, transformNumElasticIterations);

        //# Convert back to raw data
        Eigen::Map<FeatureMat>(floatingFeaturesArray, floatingFeatures.rows(), floatingFeatures.cols()) = floatingFeatures;
    }


    void downsample_mesh_mex(const float featuresArray[], const size_t numElements,
                            const int facesArray[], const size_t numFaces,
                            const float flagsArray[],
                            float sampledFeaturesArray[], const size_t numSampledElements,
                            int sampledFacesArray[], const size_t numSampledFaces,
                            float sampledFlagsArray[],
                            int originalIndicesArray[],
                            const float downsampleRatio/* = 0.8f*/){
        //# Convert arrays to Eigen matrices (see http://dovgalecs.com/blog/eigen-how-to-get-in-and-out-data-from-eigen-matrix/)
        const FeatureMat features = Eigen::Map<const FeatureMat>(featuresArray, numElements, registration::NUM_FEATURES);
        const FacesMat faces = Eigen::Map<const FacesMat>(facesArray, numFaces, 3);
        const VecDynFloat flags = Eigen::Map<const VecDynFloat>(flagsArray, numElements);
        FeatureMat sampledFeatures = Eigen::Map<FeatureMat>(sampledFeaturesArray, numSampledElements, registration::NUM_FEATURES);
        FacesMat sampledFaces = Eigen::Map<FacesMat>(sampledFacesArray, numSampledFaces, 3);
        VecDynFloat sampledFlags = Eigen::Map<VecDynFloat>(sampledFlagsArray, numSampledElements);
        VecDynInt originalIndices = Eigen::Map<VecDynInt>(originalIndicesArray, numSampledElements);

        //# Downsample
        downsample_mesh(features, faces, flags,
                        sampledFeatures, sampledFaces, sampledFlags,
                        originalIndices, downsampleRatio);

        //# Convert back to raw data
        Eigen::Map<FeatureMat>(sampledFeaturesArray, sampledFeatures.rows(), sampledFeatures.cols()) = sampledFeatures;
        Eigen::Map<FacesMat>(sampledFacesArray, sampledFaces.rows(), sampledFaces.cols()) = sampledFaces;
        Eigen::Map<VecDynFloat>(sampledFlagsArray, sampledFlags.rows(), sampledFlags.cols()) = sampledFlags;
        Eigen::Map<VecDynInt>(originalIndicesArray, originalIndices.rows(), originalIndices.cols()) = originalIndices;
    }


    void scaleshift_mesh_mex(const float oldFeaturesArray[], const size_t numOldElements,
                            const int oldIndicesArray[],
                            float newFeaturesArray[], const size_t numNewElements,
                            const int newIndicesArray[]){
        //# Convert arrays to Eigen matrices (see http://dovgalecs.com/blog/eigen-how-to-get-in-and-out-data-from-eigen-matrix/)
        const FeatureMat oldFeatures = Eigen::Map<const FeatureMat>(oldFeaturesArray, numOldElements, registration::NUM_FEATURES);
        const VecDynInt oldIndices = Eigen::Map<const VecDynInt>(oldIndicesArray, numOldElements);
        FeatureMat newFeatures = Eigen::Map<FeatureMat>(newFeaturesArray, numNewElements, registration::NUM_FEATURES);
        const VecDynInt newIndices = Eigen::Map<const VecDynInt>(newIndicesArray, numNewElements);

        //# ScaleShift
        scale_shift_mesh(oldFeatures, oldIndices,
                        newFeatures, newIndices);

        //# Convert back to raw data
        Eigen::Map<FeatureMat>(newFeaturesArray, newFeatures.rows(), newFeatures.cols()) = newFeatures;
    }

    void compute_normals_mex(const float positionsArray[], const size_t numElements,
                            const int facesArray[], const size_t numFaces,
                            float normalsArray[]){
        //# Convert arrays to Eigen matrices (see http://dovgalecs.com/blog/eigen-how-to-get-in-and-out-data-from-eigen-matrix/)
        const Vec3Mat inPositions = Eigen::Map<const Vec3Mat>(positionsArray, numElements, 3);
        const FacesMat inFaces = Eigen::Map<const FacesMat>(facesArray, numFaces, 3);
        Vec3Mat outNormals = Eigen::Map<Vec3Mat>(normalsArray, numElements, 3);

        //# ScaleShift
        compute_normals(inPositions, inFaces, outNormals);

        //# Convert back to raw data
        Eigen::Map<Vec3Mat>(normalsArray, outNormals.rows(), outNormals.cols()) = outNormals;
    }



    //######################################################################################
    //################################  REGISTRATION  ######################################
    //######################################################################################
    /*
    Full Pyramid Nonrigid Registration
    This is the function you'll normally want to call to nonrigidly register a floating mesh to a target mesh.
    */
    void pyramid_registration(FeatureMat& floatingFeatures, const FeatureMat& targetFeatures,
                                const FacesMat& floatingFaces, const FacesMat& targetFaces,
                                const VecDynFloat& floatingFlags, const VecDynFloat& targetFlags,
                                const size_t numIterations/* = 60*/, const size_t numPyramidLayers/* = 3*/,
                                const float downsampleFloatStart/* = 90*/, const float downsampleTargetStart/* = 90*/,
                                const float downsampleFloatEnd/* = 0*/, const float downsampleTargetEnd/* = 0*/,
                                const bool correspondencesSymmetric/* = true*/, const size_t correspondencesNumNeighbours/* = 5*/,
                                const float correspondencesFlagThreshold/* = 0.9f*/,
                                const float inlierKappa/* = 4.0f*/, const bool inlierUseOrientation/* = true*/,
                                const float transformSigma/* = 3.0f*/,
                                const size_t transformNumViscousIterationsStart/* = 50*/, const size_t transformNumViscousIterationsEnd/* = 1*/,
                                const size_t transformNumElasticIterationsStart/* = 50*/, const size_t transformNumElasticIterationsEnd/* = 1*/)
    {
        registration::PyramidNonrigidRegistration registrator;
        registrator.set_input(floatingFeatures, targetFeatures,
                                floatingFaces, targetFaces,
                                floatingFlags, targetFlags);
        registrator.set_parameters(numIterations, numPyramidLayers,
                                    downsampleFloatStart, downsampleTargetStart,
                                    downsampleFloatEnd, downsampleTargetEnd,
                                    correspondencesSymmetric, correspondencesNumNeighbours,
                                    correspondencesFlagThreshold,
                                    inlierKappa, inlierUseOrientation,
                                    transformSigma,
                                    transformNumViscousIterationsStart, transformNumViscousIterationsEnd,
                                    transformNumElasticIterationsStart, transformNumElasticIterationsEnd);
        registrator.update();
    }

    /*
    Standard Nonrigid Registration
    This is the standard nonrigid registration procedure without pyramid approach, so computationally a bit slower.
    */
    void nonrigid_registration(FeatureMat& floatingFeatures, const FeatureMat& targetFeatures,
                                const FacesMat& floatingFaces, const FacesMat& targetFaces,
                                const VecDynFloat& floatingFlags, const VecDynFloat& targetFlags,
                                const size_t numIterations/* = 60*/,
                                const bool correspondencesSymmetric/* = true*/, const size_t correspondencesNumNeighbours/* = 5*/,
                                const float correspondencesFlagThreshold/* = 0.9f*/,
                                const float inlierKappa/* = 4.0f*/, const bool inlierUseOrientation/* = true*/,
                                const float transformSigma/* = 3.0f*/,
                                const size_t transformNumViscousIterationsStart/* = 50*/, const size_t transformNumViscousIterationsEnd/* = 1*/,
                                const size_t transformNumElasticIterationsStart/* = 50*/, const size_t transformNumElasticIterationsEnd/* = 1*/)
    {

        registration::NonrigidRegistration registrator;
        registrator.set_input(&floatingFeatures, &targetFeatures,
                                &floatingFaces,
                                &floatingFlags, &targetFlags);
        registrator.set_parameters(correspondencesSymmetric, correspondencesNumNeighbours,
                                    correspondencesFlagThreshold,
                                    inlierKappa, inlierUseOrientation,
                                    numIterations,
                                    transformSigma,
                                    transformNumViscousIterationsStart, transformNumViscousIterationsEnd,
                                    transformNumElasticIterationsStart, transformNumElasticIterationsEnd);
        registrator.update();
    }

    /*
    Rigid Registration
    */
    void rigid_registration(FeatureMat& floatingFeatures, const FeatureMat& targetFeatures,
                                const FacesMat& floatingFaces, const FacesMat& targetFaces,
                                const VecDynFloat& floatingFlags, const VecDynFloat& targetFlags,
                                const size_t numIterations/* = 20*/,
                                const bool correspondencesSymmetric/* = true*/, const size_t correspondencesNumNeighbours/* = 5*/,
                                const float correspondencesFlagThreshold/* = 0.9f*/,
                                const float inlierKappa/* = 4.0f*/, const bool inlierUseOrientation/*=true*/)
    {
        registration::RigidRegistration registrator;
        registrator.set_input(&floatingFeatures, &targetFeatures,
                                &floatingFlags, &targetFlags);
        registrator.set_parameters(correspondencesSymmetric, correspondencesNumNeighbours,
                                    correspondencesFlagThreshold, inlierKappa,
                                    inlierUseOrientation, numIterations);
        registrator.update();
    }




    //######################################################################################
    //############################  REGISTRATION MODULES  ##################################
    //######################################################################################

    //# Correspondences
    void compute_correspondences(const FeatureMat& floatingFeatures, const FeatureMat& targetFeatures,
                                const VecDynFloat& floatingFlags, const VecDynFloat& targetFlags,
                                FeatureMat& correspondingFeatures, VecDynFloat& correspondingFlags,
                                const bool symmetric/* = true*/, const size_t numNeighbours/* = 5*/,
                                const float correspondencesFlagThreshold/* = 0.9f*/){
        registration::BaseCorrespondenceFilter* correspondenceFilter = NULL;
        if (symmetric) {
            correspondenceFilter = new registration::SymmetricCorrespondenceFilter();
        }
        else {
            correspondenceFilter = new registration::CorrespondenceFilter();
        }
        correspondenceFilter->set_floating_input(&floatingFeatures, &floatingFlags);
        correspondenceFilter->set_target_input(&targetFeatures, &targetFlags);
        correspondenceFilter->set_output(&correspondingFeatures, &correspondingFlags);
        correspondenceFilter->set_parameters(numNeighbours, correspondencesFlagThreshold);
        correspondenceFilter->update();

        delete correspondenceFilter;
    }

    //# Inliers
    void compute_inlier_weights(const FeatureMat& floatingFeatures, const FeatureMat& correspondingFeatures,
                                const VecDynFloat& correspondingFlags, VecDynFloat& inlierWeights,
                                const float kappa/* = 4.0f*/, const bool useOrientation/* = true*/){
        registration::InlierDetector inlierDetector;
        inlierDetector.set_input(&floatingFeatures, &correspondingFeatures,
                                    &correspondingFlags);
        inlierDetector.set_output(&inlierWeights);
        inlierDetector.set_parameters(kappa, useOrientation);
        inlierDetector.update();
    }

    //# Rigid Transformation
    void compute_rigid_transformation(FeatureMat& floatingFeatures, const FeatureMat& correspondingFeatures,
                                    const VecDynFloat& inlierWeights, const bool allowScaling/* = false*/){
        registration::RigidTransformer rigidTransformer;
        rigidTransformer.set_input(&correspondingFeatures, &inlierWeights);
        rigidTransformer.set_output(&floatingFeatures);
        rigidTransformer.set_parameters(allowScaling);
        rigidTransformer.update();
    }

    //# Nonrigid Transformation
    void compute_nonrigid_transformation(FeatureMat& floatingFeatures, const FeatureMat& correspondingFeatures,
                                        const FacesMat& floatingFaces, const VecDynFloat& floatingFlags,
                                        const VecDynFloat& inlierWeights,
                                        const size_t numSmoothingNeighbours/* = 10*/, const float sigmaSmoothing/* = 3.0f*/,
                                        const size_t numViscousIterations/* = 50*/, const size_t numElasticIterations/* = 50*/){
        registration::ViscoElasticTransformer transformer;
        transformer.set_input(&correspondingFeatures, &inlierWeights, &floatingFlags, &floatingFaces);
        transformer.set_output(&floatingFeatures);
        transformer.set_parameters(numSmoothingNeighbours, sigmaSmoothing, numViscousIterations,numElasticIterations);
        transformer.update();
    }


    //# Downsampler
    void downsample_mesh(const FeatureMat& features, const FacesMat& faces,
                        const VecDynFloat& flags,
                        FeatureMat& downsampledFeatures, FacesMat& downsampledFaces,
                        VecDynFloat& downsampledFlags, VecDynInt& originalIndices,
                        const float downsampleRatio/* = 0.8f*/){

        registration::Downsampler downsampler;
        downsampler.set_input(&features, &faces, &flags);
        downsampler.set_output(downsampledFeatures, downsampledFaces,
                                downsampledFlags, originalIndices);
        downsampler.set_parameters(downsampleRatio);
        downsampler.update();
    }


    //# ScaleShifter
    //## The scaleshifter is meant to transition from one scale in the pyramid to the next.
    void scale_shift_mesh(const FeatureMat& previousFeatures, const VecDynInt& previousIndices,
                        FeatureMat& newFeatures, const VecDynInt& newIndices){
        registration::ScaleShifter scaleShifter;
        scaleShifter.set_input(previousFeatures, previousIndices, newIndices);
        scaleShifter.set_output(newFeatures);
        scaleShifter.update();
    }


    //######################################################################################
    //###############################  MESH OPERATIONS  ####################################
    //######################################################################################
    //# Compute Normals from positions and faces
    void compute_normals(const Vec3Mat &inPositions, const FacesMat &inFaces,
                        Vec3Mat &outNormals){
        registration::update_normals_for_altered_positions(inPositions, inFaces, outNormals);
    }

    //######################################################################################
    //################################  INPUT/OUTPUT  ######################################
    //######################################################################################

    void read_obj_files(const std::string floatingMeshPath, const std::string targetMeshPath,
                        FeatureMat& floatingFeatures, FeatureMat& targetFeatures,
                        FacesMat& floatingFaces, FacesMat& targetFaces){
        registration::import_data(floatingMeshPath, targetMeshPath,
                                    floatingFeatures, targetFeatures,
                                    floatingFaces, targetFaces);
    }

    void write_obj_files(FeatureMat& features, FacesMat& faces, const std::string meshPath){
        registration::export_data(features, faces, meshPath);
    }

#ifdef __cplusplus
}//extern C
#endif // __cplusplus

}//namespace meshmonk
