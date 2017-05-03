#ifndef MESHMONK_HPP
#define MESHMONK_HPP

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/reader/OBJReader.hh>
#include <OpenMesh/Core/IO/writer/OBJWriter.hh>
#include <OpenMesh/Tools/Decimater/DecimaterT.hh>
#include <OpenMesh/Tools/Decimater/ModQuadricT.hh>
#include <Eigen/Dense>
#include "src/PyramidNonrigidRegistration.hpp"
#include "src/RigidRegistration.hpp"
#include "src/NonrigidRegistration.hpp"
#include "src/InlierDetector.hpp"
#include "src/CorrespondenceFilter.hpp"
#include "src/SymmetricCorrespondenceFilter.hpp"
#include "src/RigidTransformer.hpp"
#include "src/ViscoElasticTransformer.hpp"
#include "src/Downsampler.hpp"
#include "src/ScaleShifter.hpp"
#include "global.hpp"
#include "src/helper_functions.hpp"



typedef OpenMesh::DefaultTraits MyTraits;
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>  TriMesh;
typedef OpenMesh::Decimater::DecimaterT<TriMesh>    DecimaterType;
typedef OpenMesh::Decimater::ModQuadricT<TriMesh>::Handle HModQuadric;
typedef Eigen::Matrix< int, Eigen::Dynamic, 3> FacesMat; //matrix Mx3 of type unsigned int
typedef Eigen::VectorXf VecDynFloat;
typedef Eigen::Matrix< float, Eigen::Dynamic, registration::NUM_FEATURES> FeatureMat; //matrix Mx6 of type float
typedef Eigen::MatrixX3f Vec3Mat;


namespace meshmonk{

#ifdef __cplusplus
extern "C"
#endif // __cplusplus
{
    //######################################################################################
    //################################  TEST SHIZZLE  ######################################
    //######################################################################################
    /*
    We're implementing this function simply to test MEX'ing in MATLAB.
    */
    void test_meshmonk_mexing(FeatureMat& floatingFeatures, const FeatureMat& targetFeatures, const float multiplier = 2.0f){
        floatingFeatures += targetFeatures * multiplier;
    }

    /*
    Raw data version of test_meshmonk_mexing()
    */
    void test_meshmonk_mexing_raw(float floatingFeaturesRaw[], const float targetFeaturesRaw[],
                                    const size_t numFloatingElements, const size_t numTargetElements,
                                    const float multiplier = 2.0f){
        //# Convert raw data pointers to Eigen matrices (see http://dovgalecs.com/blog/eigen-how-to-get-in-and-out-data-from-eigen-matrix/)
        FeatureMat floatingFeatures = Eigen::Map<FeatureMat>(floatingFeaturesRaw, numFloatingElements, registration::NUM_FEATURES);
        const FeatureMat targetFeatures = Eigen::Map<const FeatureMat>(targetFeaturesRaw, numTargetElements, registration::NUM_FEATURES);

        //# Call test_meshmonk_mexing()
        test_meshmonk_mexing(floatingFeatures, targetFeatures, multiplier);
        floatingFeatures(0,5) *= 2.0f;
        floatingFeatures(5,1) *= 3.0f;

        //# Convert back to raw data
        Eigen::Map<FeatureMat>(floatingFeaturesRaw, floatingFeatures.rows(), floatingFeatures.cols()) = floatingFeatures;
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
                                const size_t numIterations = 60, const size_t numPyramidLayers = 3,
                                const float downsampleFloatStart = 90, const float downsampleTargetStart = 90,
                                const float downsampleFloatEnd = 0, const float downsampleTargetEnd = 0,
                                const bool correspondencesSymmetric = true, const size_t correspondencesNumNeighbours = 5,
                                const float inlierKappa = 4.0f, const bool inlierUseOrientation = true,
                                const float transformSigma = 3.0f,
                                const size_t transformNumViscousIterationsStart = 50, const size_t transformNumViscousIterationsEnd = 1,
                                const size_t transformNumElasticIterationsStart = 50, const size_t transformNumElasticIterationsEnd = 1);

    /*
    Standard Nonrigid Registration
    This is the standard nonrigid registration procedure without pyramid approach, so computationally a bit slower.
    */
    void nonrigid_registration(FeatureMat& floatingFeatures, const FeatureMat& targetFeatures,
                                const FacesMat& floatingFaces, const FacesMat& targetFaces,
                                const VecDynFloat& floatingFlags, const VecDynFloat& targetFlags,
                                const size_t numIterations = 60,
                                const bool correspondencesSymmetric = true, const size_t correspondencesNumNeighbours = 5,
                                const float inlierKappa = 4.0f, const bool inlierUseOrientation = true,
                                const float transformSigma = 3.0f,
                                const size_t transformNumViscousIterationsStart = 50, const size_t transformNumViscousIterationsEnd = 1,
                                const size_t transformNumElasticIterationsStart = 50, const size_t transformNumElasticIterationsEnd = 1);

    /*
    Rigid Registration
    */
    void rigid_registration(FeatureMat& floatingFeatures, const FeatureMat& targetFeatures,
                                const FacesMat& floatingFaces, const FacesMat& targetFaces,
                                const VecDynFloat& floatingFlags, const VecDynFloat& targetFlags,
                                const size_t numIterations = 20,
                                const bool correspondencesSymmetric = true, const size_t correspondencesNumNeighbours = 5,
                                const float inlierKappa = 4.0f, const bool inlierUseOrientation = true);




    //######################################################################################
    //############################  REGISTRATION MODULES  ##################################
    //######################################################################################

    //# Correspondences
    void compute_correspondences(const FeatureMat& floatingFeatures, const FeatureMat& targetFeatures,
                                const VecDynFloat& floatingFlags, const VecDynFloat& targetFlags,
                                FeatureMat& correspondingFeatures, VecDynFloat& correspondingFlags,
                                const bool symmetric = true, const size_t numNeighbours = 5);

    //# Inliers
    void compute_inlier_weights(const FeatureMat& floatingFeatures, const FeatureMat& correspondingFeatures,
                                const VecDynFloat& correspondingFlags, VecDynFloat& inlierWeights,
                                const float kappa = 4.0f, const bool useOrientation = true);

    //# Rigid Transformation
    void compute_rigid_transformation(FeatureMat& floatingFeatures, const FeatureMat& correspondingFeatures,
                                    const VecDynFloat& inlierWeights, const bool allowScaling = false);

    //# Nonrigid Transformation
    void compute_nonrigid_transformation(FeatureMat& floatingFeatures, const FeatureMat& correspondingFeatures,
                                        const FacesMat& floatingFaces, const VecDynFloat& floatingFlags,
                                        const VecDynFloat& inlierWeights,
                                        const size_t numSmoothingNeighbours = 10, const float sigmaSmoothing = 3.0f,
                                        const size_t numViscousIterations = 50, const size_t numElasticIterations = 50);


    //# Downsampler
    void downsample_mesh(const FeatureMat& features, const FacesMat& faces,
                        const VecDynFloat& flags,
                        FeatureMat& downsampledFeatures, FacesMat& downsampledFaces,
                        VecDynFloat& downsampledFlags, VecDynInt& originalIndices,
                        const float downsampleRatio = 0.8f);


    //# ScaleShifter
    //## The scaleshifter is meant to transition from one scale in the pyramid to the next.
    void scale_shift_mesh(const FeatureMat& previousFeatures, const VecDynInt& previousIndices,
                        FeatureMat& newFeatures, const VecDynInt& newIndices);

    //######################################################################################
    //###############################  MESH OPERATIONS  ####################################
    //######################################################################################
    void compute_normals(const Vec3Mat &inPositions, const FacesMat &inFaces,
                        Vec3Mat &outNormals);


    //######################################################################################
    //################################  INPUT/OUTPUT  ######################################
    //######################################################################################

    void read_obj_files(const std::string floatingMeshPath, const std::string targetMeshPath,
                        FeatureMat& floatingFeatures, FeatureMat& targetFeatures,
                        FacesMat& floatingFaces, FacesMat& targetFaces);

    void write_obj_files(FeatureMat& features, FacesMat& faces, const std::string meshPath);


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
                                const size_t numIterations = 60, const size_t numPyramidLayers = 3,
                                const float downsampleFloatStart = 90, const float downsampleTargetStart = 90,
                                const float downsampleFloatEnd = 0, const float downsampleTargetEnd = 0,
                                const bool correspondencesSymmetric = true, const size_t correspondencesNumNeighbours = 5,
                                const float inlierKappa = 4.0f, const bool inlierUseOrientation = true,
                                const float transformSigma = 3.0f,
                                const size_t transformNumViscousIterationsStart = 50, const size_t transformNumViscousIterationsEnd = 1,
                                const size_t transformNumElasticIterationsStart = 50, const size_t transformNumElasticIterationsEnd = 1);

    void nonrigid_registration_mex(float floatingFeaturesArray[], const float targetFeaturesArray[],
                                const size_t numFloatingElements, const size_t numTargetElements,
                                const int floatingFacesArray[], const int targetFacesArray[],
                                const size_t numFloatingFaces, const size_t numTargetFaces,
                                const float floatingFlagsArray[], const float targetFlagsArray[],
                                const size_t numIterations = 60,
                                const bool correspondencesSymmetric = true, const size_t correspondencesNumNeighbours = 5,
                                const float inlierKappa = 4.0f, const bool inlierUseOrientation = true,
                                const float transformSigma = 3.0f,
                                const size_t transformNumViscousIterationsStart = 50, const size_t transformNumViscousIterationsEnd = 1,
                                const size_t transformNumElasticIterationsStart = 50, const size_t transformNumElasticIterationsEnd = 1);

    void rigid_registration_mex(float floatingFeaturesArray[], const float targetFeaturesArray[],
                                const size_t numFloatingElements, const size_t numTargetElements,
                                const int floatingFacesArray[], const int targetFacesArray[],
                                const size_t numFloatingFaces, const size_t numTargetFaces,
                                const float floatingFlagsArray[], const float targetFlagsArray[],
                                const size_t numIterations = 60,
                                const bool correspondencesSymmetric = true, const size_t correspondencesNumNeighbours = 5,
                                const float inlierKappa = 4.0f, const bool inlierUseOrientation = true);

    void compute_correspondences_mex(const float floatingFeaturesArray[], const float targetFeaturesArray[],
                                    const size_t numFloatingElements, const size_t numTargetElements,
                                    const float floatingFlagsArray[], const float targetFlagsArray[],
                                    float correspondingFeaturesArray[], float correspondingFlagsArray[],
                                    const bool correspondencesSymmetric = true, const size_t correspondencesNumNeighbours = 5);

    void compute_inlier_weights_mex(const float floatingFeaturesArray[], const float correspondingFeaturesArray[],
                                    const size_t numFloatingElements,
                                    const float correspondingFlagsArray[], float inlierWeightsArray[],
                                    const float inlierKappa/*= 4.0f*/, const bool useOrientation/*= true*/);

    void compute_rigid_transformation_mex(float floatingFeaturesArray[], const size_t numFloatingElements,
                                        const float correspondingFeaturesArray[], const float inlierWeightsArray[],
                                        const bool allowScaling /*= fakse*/);

    void compute_nonrigid_transformation_mex(float floatingFeaturesArray[], const float correspondingFeaturesArray[],
                                            const size_t numFloatingElements,
                                            const int floatingFacesArray[], const size_t numFloatingFaces,
                                            const float floatingFlagsArray[], const float inlierWeightsArray[],
                                            const size_t transformNumNeighbours/*= 10*/, const float transformSigma/*= 3.0f*/,
                                            const size_t transformNumViscousIterations/*= 50*/, const size_t transformNumElasticIterations/*= 50*/);

    void downsample_mesh_mex(const float featuresArray[], const size_t numElements,
                            const int facesArray[], const size_t numFaces,
                            const float flagsArray[],
                            float sampledFeaturesArray[], const size_t numSampledElements,
                            int sampledFacesArray[], const size_t numSampledFaces,
                            float sampledFlagsArray[],
                            int originalIndicesArray[],
                            const float downsampleRatio/* = 0.8f*/);

    void scaleshift_mesh_mex(const float oldFeaturesArray[], const size_t numOldElements,
                            const int oldIndicesArray[],
                            float newFeaturesArray[], const size_t numNewElements,
                            const int newIndicesArray[]);

    void compute_normals_mex(const float positionsArray[], const size_t numElements,
                            const int facesArray[], const size_t numFaces,
                            float normalsArray[]);

#ifdef __cplusplus
}//extern C
#endif // __cplusplus

}//namespace meshmonk

#endif //MESHMONK_HPP
