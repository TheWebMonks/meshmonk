#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/reader/OBJReader.hh>
#include <OpenMesh/Core/IO/writer/OBJWriter.hh>
#include <iostream>
#include <Eigen/Dense>
#include <stdio.h>
#include <math.h>
#include "meshmonk.hpp"

using namespace std;

typedef OpenMesh::DefaultTraits MyTraits;
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>  TriMesh;
typedef Eigen::Matrix< int, Eigen::Dynamic, 3> FacesMat; //matrix Mx3 of type unsigned int
typedef Eigen::VectorXf VecDynFloat;
typedef Eigen::Matrix< float, Eigen::Dynamic, 6> FeatureMat; //matrix Mx6 of type float


int main()
{

    /*
    ############################################################################
    ##############################  INPUT  #####################################
    ############################################################################
    */
    //# IO variables
    //const string floatingDir = "/home/jonatan/projects/meshmonk/examples/data/karlijne/Normaal/Cleaned/CL0009A_Clean.obj";
//    const string floatingDir = "/home/jonatan/projects/meshmonk/examples/data/karlijne/Normaal/Cleaned/rigidResult.obj";
//    const string targetDir = "/home/jonatan/projects/meshmonk/examples/data/karlijne/Normaal/Cleaned/CL0013A_Clean.obj";
//    const string floatingDir = "/home/jonatan/projects/meshmonk/examples/faceTemplate.obj";
//    const string targetDir = "/home/jonatan/projects/meshmonk/examples/faceTarget.obj";
//    const string resultDir = "/home/jonatan/projects/meshmonk/examples/data/karlijne/Normaal/Cleaned/resultexample.obj";
//    const string floatingDir = "/home/jonatan/projects/meshmonk/examples/data/bunny.obj";
//    const string targetDir = "/home/jonatan/projects/meshmonk/examples/data/bunny2.obj";
//    const string resultDir = "/home/jonatan/projects/meshmonk/examples/data/bunnyResult.obj";

    //# My Matlab Data
//    const string floatingDir = "/home/jonatan/projects/meshmonk/examples/faceTemplate.obj";
//    const string targetDir = "/home/jonatan/projects/meshmonk/examples/faceTarget.obj";

    //# Peter's Matlab Data
    const string floatingDir = "/Users/maartenhuijsmans/GitProjects/meshmonk/meshmonk/floating.obj";
    const string targetDir = "/Users/maartenhuijsmans/GitProjects/meshmonk/meshmonk/target.obj";



    const string resultDir = "/home/jonatan/projects/meshmonk/examples/matlabResult.obj";





    //# Load meshes and convert to feature matrices
    FeatureMat floatingFeatures;
    FeatureMat targetFeatures;
    FacesMat floatingFaces;
    FacesMat targetFaces;
    meshmonk::read_obj_files(floatingDir, targetDir,
                    floatingFeatures, targetFeatures,
                    floatingFaces, targetFaces);






    const size_t numFloatingVertices = floatingFeatures.rows();
    const size_t numTargetVertices = targetFeatures.rows();
    VecDynFloat floatingFlags = VecDynFloat::Ones(numFloatingVertices);
    VecDynFloat targetFlags = VecDynFloat::Ones(numTargetVertices);
//
//
//
//    /*
//    ############################################################################
//    ##############################  REGISTRATION  ##############################
//    ############################################################################
////    */
    //# Parameters
    //## Pyramid
    const size_t numIterations = 60;
    const size_t numRigidIterations = 20;
    const size_t numPyramidLayers = 1;
    const float downsampleFloatStart = 0;
    const float downsampleTargetStart = 0;
    const float downsampleFloatEnd = 0;
    const float downsampleTargetEnd = 0;

    //## Correspondences
    const bool correspondencesSymmetric = true;
    const size_t correspondencesNumNeighbours = 5;

    //## Inliers
    const float inlierKappa = 4.0f;
    const bool inlierUseOrientation = false;

    //## Transformation
    const float transformSigma = 3.0f;
    const size_t transformNumViscousIterationsStart = 55;
    const size_t transformNumViscousIterationsEnd = 2;
    const size_t transformNumElasticIterationsStart = 56;
    const size_t transformNumElasticIterationsEnd = 3;

    //# Rigid registration
    //meshmonk::rigid_registration(floatingFeatures, targetFeatures, floatingFaces, targetFaces,
//                        floatingFlags, targetFlags, numRigidIterations,
//                        correspondencesSymmetric, correspondencesNumNeighbours, inlierKappa);
//
//
    //DEBUG
        std::cout << " -- Floating Features --\n" << floatingFeatures.topRows(10) << std::endl;
        std::cout << " -- Target Features --\n" << targetFeatures.topRows(10) << std::endl;
        std::cout << " -- Floating Faces --\n" << floatingFaces.topRows(10) << std::endl;
        std::cout << " -- Target Faces --\n" << targetFaces.topRows(10) << std::endl;
        std::cout << " -- Floating Flags --\n" << floatingFlags.topRows(10) << std::endl;
        std::cout << " -- Target Flags --\n" << targetFlags.topRows(10) << std::endl;
        //END DEBUG

        //# Call pyramid_registration()
        //DEBUG
        std::cout << "Num Floating Elements              - " << floatingFeatures.rows() << std::endl;
        std::cout << "Num Target Elements                - " << targetFeatures.rows() << std::endl;
        std::cout << "Num Floating Faces                 - " << floatingFaces.rows() << std::endl;
        std::cout << "Num Target Faces                   - " << targetFaces.rows() << std::endl;
        std::cout << "Num Iterations                     - " << numIterations << std::endl;
        std::cout << "Num Pyramid Layers                 - " << numPyramidLayers << std::endl;
        std::cout << "Downsample Float Start             - " << downsampleFloatStart << std::endl;
        std::cout << "Downsample Target Start            - " << downsampleTargetStart << std::endl;
        std::cout << "Downsample Float End               - " << downsampleFloatEnd << std::endl;
        std::cout << "Downsample Target End              - " << downsampleTargetEnd << std::endl;
        std::cout << "Symmetric Correspondences          - " << correspondencesSymmetric << std::endl;
        std::cout << "Num Neighbours                     - " << correspondencesNumNeighbours << std::endl;
        std::cout << "Inlier Kappa                       - " << inlierKappa << std::endl;
        std::cout << "Transform Sigma                    - " << transformSigma << std::endl;
        std::cout << "transformNumViscousIterationsStart - " << transformNumViscousIterationsStart << std::endl;
        std::cout << "transformNumViscousIterationsEnd   - " << transformNumViscousIterationsEnd << std::endl;
        std::cout << "transformNumElasticIterationsStart - " << transformNumElasticIterationsStart << std::endl;
        std::cout << "transformNumElasticIterationsEnd   - " << transformNumElasticIterationsEnd << std::endl;
        //END DEBUG
//
//
    //# Nonrigid registration
    //DEBUG
    //## Initialize variables
//    floatingFeatures = FeatureMat::Zero(3,6);
//    floatingFeatures << 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,
//                        1.0, 1.0, 0.0, 1.0, 0.0, 0.0,
//                        1.0, 0.0, 0.0, 1.0, 0.0, 0.0;
//    targetFeatures = FeatureMat::Zero(3,6);
//    targetFeatures << 1.1, 1.0, 1.0, 1.0, 0.0, 0.0,
//                      0.9, 1.0, 0.0, 1.0, 0.0, 0.0,
//                      1.2, 0.0, 0.0, 1.0, 0.0, 0.0;
//    floatingFaces = FacesMat::Zero(1,3);
//    floatingFaces << 0, 1, 2;
//    targetFaces = FacesMat::Zero(1,3);
//    targetFaces << 0, 1, 2;
//    floatingFlags = VecDynFloat::Ones(3);
//    targetFlags = VecDynFloat::Ones(3);
//
//
//    std::cout << "floatingFeatures: " << floatingFeatures << std::endl;
//    std::cout << "targetFeatures: " << targetFeatures << std::endl;
//    std::cout << "floatingFaces: " << floatingFaces << std::endl;
//    std::cout << "targetFaces: " << targetFaces << std::endl;
//    std::cout << "floatingFlags: " << floatingFlags << std::endl;
//    std::cout << "targetFlags: " << targetFlags << std::endl;
//
//    //## Convert to arrays
//    float floatingFeaturesArray[18];
//    Eigen::Map<FeatureMat>(floatingFeaturesArray, floatingFeatures.rows(), floatingFeatures.cols()) = floatingFeatures;
//    float targetFeaturesArray[18];
//    Eigen::Map<FeatureMat>(targetFeaturesArray, targetFeatures.rows(), targetFeatures.cols()) = targetFeatures;
//
//    int floatingFacesArray[18];
//    Eigen::Map<FacesMat>(floatingFacesArray, floatingFaces.rows(), floatingFaces.cols()) = floatingFaces;
//    int targetFacesArray[18];
//    Eigen::Map<FacesMat>(targetFacesArray, targetFaces.rows(), targetFaces.cols()) = targetFaces;
//
//    float floatingFlagsArray[18];
//    Eigen::Map<VecDynFloat>(floatingFlagsArray, floatingFlags.rows(), floatingFlags.cols()) = floatingFlags;
//    float targetFlagsArray[18];
//    Eigen::Map<VecDynFloat>(targetFlagsArray, targetFlags.rows(), targetFlags.cols()) = targetFlags;
//
//
//
//    std::cout << "floatingFeaturesArray: " << std::endl;
//    for (int i = 0 ; i < 3 ; i++){
//        for (int j = 0 ; j < 6 ; j++) {
//            std::cout << floatingFeaturesArray[6*i + j] << " | ";
//        }
//        std::cout << std::endl;
//    }
//
//    std::cout << "targetFeaturesArray: " << targetFeaturesArray << std::endl;
//    std::cout << "floatingFacesArray: " << floatingFacesArray << std::endl;
//    std::cout << "targetFacesArray: " << targetFacesArray << std::endl;
//    std::cout << "floatingFlagsArray: " << floatingFlagsArray << std::endl;
//    std::cout << "targetFlagsArray: " << targetFlagsArray << std::endl;
    //END DEBUG

    meshmonk::pyramid_registration(floatingFeatures, targetFeatures,
                                floatingFaces, targetFaces,
                                floatingFlags, targetFlags,
                                numIterations, numPyramidLayers,
                                downsampleFloatStart, downsampleTargetStart,
                                downsampleFloatEnd, downsampleTargetEnd,
                                correspondencesSymmetric, correspondencesNumNeighbours,
                                inlierKappa, inlierUseOrientation,
                                transformSigma,
                                transformNumViscousIterationsStart, transformNumViscousIterationsEnd,
                                transformNumElasticIterationsStart, transformNumElasticIterationsEnd);
//    meshmonk::nonrigid_registration_mex(floatingFeaturesArray, targetFeaturesArray,
//                                3, 3,
//                                floatingFacesArray, targetFacesArray,
//                                1, 1,
//                                floatingFlagsArray, targetFlagsArray,
//                                numIterations,
//                                correspondencesSymmetric, correspondencesNumNeighbours,
//                                inlierKappa, inlierUseOrientation,
//                                transformSigma,
//                                transformNumViscousIterationsStart, transformNumViscousIterationsEnd,
//                                transformNumElasticIterationsStart, transformNumElasticIterationsEnd);

//    /*
//    ############################################################################
//    ##############################  OUTPUT #####################################
//    ############################################################################
//    */
//    //# Write result to file
//    std::cout << "Writing results to file..." << std::endl;
//    //DEBUG
//    std::cout << "Num Floating Elements              - " << floatingFeatures.rows() << std::endl;
//    std::cout << "Num Target Elements                - " << targetFeatures.rows() << std::endl;
//    std::cout << "Num Floating Faces                 - " << floatingFaces.rows() << std::endl;
//    std::cout << "Num Target Faces                   - " << targetFaces.rows() << std::endl;
//    //END DEBUG
//    meshmonk::write_obj_files(floatingFeatures,floatingFaces, resultDir);
//    std::cout << "Process finished." << std::endl;


    return 0;
}
