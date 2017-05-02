#include "mex.h"
#include "/home/jonatan/projects/meshmonk/meshmonk/meshmonk.hpp"
#include "mystream.cpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    //# Check input
    //## Number of input arguments
    if(nlhs != 0) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs",
                      "Zero LHS output required.");
    }
    //## Number of output arguments
    if(nrhs != 15) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                      "20 inputs required.");
    }
    
    //# Get Inputs
    //## Floating Features
    float *floatingFeatures = reinterpret_cast<float *>(mxGetData(prhs[0]));
    mwSize numFloatingElements = mxGetM(prhs[0]);
    std::cout << "Num Floating Elements - " << numFloatingElements << std::endl;
    //## Target Features
    float *targetFeatures = reinterpret_cast<float *>(mxGetData(prhs[1]));
    mwSize numTargetElements = mxGetM(prhs[1]);
    std::cout << "Num Target Elements   - " << numTargetElements << std::endl;
    //## Floating Faces
    int *floatingFaces = reinterpret_cast<int *>(mxGetData(prhs[2]));
    mwSize numFloatingFaces = mxGetM(prhs[2]);
    std::cout << "Num Floating Faces    - " << numFloatingFaces << std::endl;
    //## Target Faces
    int *targetFaces = reinterpret_cast<int *>(mxGetData(prhs[3]));
    mwSize numTargetFaces = mxGetM(prhs[3]);
    std::cout << "Num Target Faces      - " << numTargetFaces << std::endl;
    //## FLoating Flags
    float *floatingFlags = reinterpret_cast<float *>(mxGetData(prhs[4]));
    //## Target Flags
    float *targetFlags = reinterpret_cast<float *>(mxGetData(prhs[5]));
    //## Parameters
    //### Total number of iterations
    mwSize numIterations = static_cast<mwSize>(mxGetScalar(prhs[6]));
    std::cout << "Num Iterations                     - " << numIterations << std::endl;
    //### Use symmetric correspondences
    bool correspondencesSymmetric = static_cast<bool>(mxGetScalar(prhs[7]));
    std::cout << "Symmetric Correspondences          - " << correspondencesSymmetric << std::endl;
    //### Number of neighbours to use to compute corresponding points
    mwSize correspondencesNumNeighbours = static_cast<mwSize>(mxGetScalar(prhs[8]));
    std::cout << "Num Neighbours                     - " << correspondencesNumNeighbours << std::endl;
    //### Inlier kappa
    float inlierKappa = static_cast<float>(mxGetScalar(prhs[9]));
    std::cout << "Inlier Kappa                       - " << inlierKappa << std::endl;
    //### Sigma of gaussian used in visco-elastic smoothing of deformation field
    float transformSigma = static_cast<float>(mxGetScalar(prhs[10]));
    std::cout << "Transform Sigma                    - " << transformSigma << std::endl;
    //### Starting number of viscous smoothing iterations
    mwSize transformNumViscousIterationsStart = static_cast<mwSize>(mxGetScalar(prhs[11]));
    std::cout << "transformNumViscousIterationsStart - " << transformNumViscousIterationsStart << std::endl;
    //### Final number of viscous smoothing iterations
    mwSize transformNumViscousIterationsEnd = static_cast<mwSize>(mxGetScalar(prhs[12]));
    std::cout << "transformNumViscousIterationsEnd   - " << transformNumViscousIterationsEnd << std::endl;
    //### Starting number of elastic smoothing iterations
    mwSize transformNumElasticIterationsStart = static_cast<mwSize>(mxGetScalar(prhs[13]));
    std::cout << "transformNumElasticIterationsStart - " << transformNumElasticIterationsStart << std::endl;
    //### Final number of elastic smoothing iterations
    mwSize transformNumElasticIterationsEnd = static_cast<mwSize>(mxGetScalar(prhs[14]));
    std::cout << "transformNumElasticIterationsEnd   - " << transformNumElasticIterationsEnd << std::endl;
    
    //# Execute c++ function
//     meshmonk::test_meshmonk_mexing_raw(floatingFeatures, targetFeatures, numFloatingElements, numTargetElements, transformSigma);
    meshmonk::nonrigid_registration_mex(floatingFeatures, targetFeatures,
                                numFloatingElements, numTargetElements,
                                floatingFaces, targetFaces,
                                numFloatingFaces, numTargetFaces,
                                floatingFlags, targetFlags,
                                numIterations,
                                correspondencesSymmetric, correspondencesNumNeighbours,
                                inlierKappa,
                                transformSigma,
                                transformNumViscousIterationsStart, transformNumViscousIterationsEnd,
                                transformNumElasticIterationsStart, transformNumElasticIterationsEnd);
    
//     //# Set Output
//     int numCols = 6;
//     plhs[0] = mxCreateNumericMatrix(numFloatingElements, numCols, mxSINGLE_CLASS, mxREAL); // output: double matrix
//     auto output = mxGetPr(plhs[0]);
//     //## Copy result form c++ function into the output
//     for (unsigned i = 0 ; i < numFloatingElements * numCols ; i++){
//         output[i] = floatingFeatures[i];
//     }
  
}