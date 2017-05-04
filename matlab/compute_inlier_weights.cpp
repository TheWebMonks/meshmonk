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
    if(nrhs != 6) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                      "6 inputs required.");
    }
    
    //# Get Inputs
    //## Floating Features
    float *floatingFeatures = reinterpret_cast<float *>(mxGetData(prhs[0]));
    mwSize numFloatingElements = mxGetM(prhs[0]);
    std::cout << "Num Floating Elements - " << numFloatingElements << std::endl;
    //## Corresponding Features
    float *correspondingFeatures = reinterpret_cast<float *>(mxGetData(prhs[1]));
    //## Corresponding Flags
    float *correspondingFlags = reinterpret_cast<float *>(mxGetData(prhs[2]));
    //## Inlier Weights
    float *inlierWeights = reinterpret_cast<float *>(mxGetData(prhs[3]));
    //## Parameters
    //### Inlier kappa
    float inlierKappa = static_cast<float>(mxGetScalar(prhs[4]));
    std::cout << "Inlier Kappa                       - " << inlierKappa << std::endl;
    //### Inlier Orientation
    float inlierUseOrientation = static_cast<float>(mxGetScalar(prhs[5]));
    std::cout << "Inlier Use Orientation             - " << inlierUseOrientation << std::endl;
    
    //# Execute c++ function
    meshmonk::compute_inlier_weights_mex(floatingFeatures, correspondingFeatures,
                                        numFloatingElements,
                                        correspondingFlags, inlierWeights,
                                        inlierKappa, inlierUseOrientation);
    
//     //# Set Output
//     int numCols = 1;
//     plhs[0] = mxCreateNumericMatrix(numFloatingElements, numCols, mxSINGLE_CLASS, mxREAL);
//     auto output = mxGetPr(plhs[0]);
//     //## Copy result form c++ function into the output
//     for (unsigned i = 0 ; i < numFloatingElements * numCols ; i++){
//         output[i] = inlierWeights[i];
//     }
  
}