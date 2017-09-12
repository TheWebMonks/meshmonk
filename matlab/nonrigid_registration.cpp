#include "mex.h"
#include <meshmonk/meshmonk.hpp>
#include "mystream.cpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    //# Check input
    //## Number of input arguments
    if(nlhs != 0) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs",
                      "Zero LHS output required.");
    }
    //## Number of output arguments
    if(nrhs != 18) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                      "18 inputs required.");
    }
    
    //# Get Inputs
    //## Floating Features
    float *floatingFeatures = reinterpret_cast<float *>(mxGetData(prhs[0]));
    mwSize numFloatingElements = mxGetM(prhs[0]);
    //## Target Features
    float *targetFeatures = reinterpret_cast<float *>(mxGetData(prhs[1]));
    mwSize numTargetElements = mxGetM(prhs[1]);
    //## Floating Faces
    int *floatingFaces = reinterpret_cast<int *>(mxGetData(prhs[2]));
    mwSize numFloatingFaces = mxGetM(prhs[2]);
    //## Target Faces
    int *targetFaces = reinterpret_cast<int *>(mxGetData(prhs[3]));
    mwSize numTargetFaces = mxGetM(prhs[3]);
    //## FLoating Flags
    float *floatingFlags = reinterpret_cast<float *>(mxGetData(prhs[4]));
    //## Target Flags
    float *targetFlags = reinterpret_cast<float *>(mxGetData(prhs[5]));
    //## Parameters
    //### Total number of iterations
    mwSize numIterations = static_cast<mwSize>(mxGetScalar(prhs[6]));
    //### Use symmetric correspondences
    bool correspondencesSymmetric = static_cast<bool>(mxGetScalar(prhs[7]));
    //### Number of neighbours to use to compute corresponding points
    mwSize correspondencesNumNeighbours = static_cast<mwSize>(mxGetScalar(prhs[8]));
    //### Flag threshold to mark corresponding flag as 0.0 or 1.0
    float correspondencesFlagThreshold = static_cast<float>(mxGetScalar(prhs[9]));
    //### Equalize the push and pull forces (when using symmetric correspondences)
    bool correspondencesEqualizePushPull = static_cast<bool>(mxGetScalar(prhs[10]));
    //### Inlier kappa
    float inlierKappa = static_cast<float>(mxGetScalar(prhs[11]));
    //### Inlier Orientation
    float inlierUseOrientation = static_cast<float>(mxGetScalar(prhs[12]));
    //### Sigma of gaussian used in visco-elastic smoothing of deformation field
    float transformSigma = static_cast<float>(mxGetScalar(prhs[13]));
    //### Starting number of viscous smoothing iterations
    mwSize transformNumViscousIterationsStart = static_cast<mwSize>(mxGetScalar(prhs[14]));
    //### Final number of viscous smoothing iterations
    mwSize transformNumViscousIterationsEnd = static_cast<mwSize>(mxGetScalar(prhs[15]));
    //### Starting number of elastic smoothing iterations
    mwSize transformNumElasticIterationsStart = static_cast<mwSize>(mxGetScalar(prhs[16]));
    //### Final number of elastic smoothing iterations
    mwSize transformNumElasticIterationsEnd = static_cast<mwSize>(mxGetScalar(prhs[17]));
    
    //# Execute c++ function
//     meshmonk::test_meshmonk_mexing_raw(floatingFeatures, targetFeatures, numFloatingElements, numTargetElements, transformSigma);
    meshmonk::nonrigid_registration_mex(floatingFeatures, targetFeatures,
                                numFloatingElements, numTargetElements,
                                floatingFaces, targetFaces,
                                numFloatingFaces, numTargetFaces,
                                floatingFlags, targetFlags,
                                numIterations,
                                correspondencesSymmetric, correspondencesNumNeighbours,
                                correspondencesFlagThreshold, correspondencesEqualizePushPull,
                                inlierKappa, inlierUseOrientation,
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