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
    if(nrhs != 5) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                      "5 inputs required.");
    }
    
    //# Get Inputs
    //## Floating Features
    float *floatingFeatures = reinterpret_cast<float *>(mxGetData(prhs[0]));
    mwSize numFloatingElements = mxGetM(prhs[0]);
    //## Corresponding Features
    float *correspondingFeatures = reinterpret_cast<float *>(mxGetData(prhs[1]));
    //## Inlier Weights
    float *inlierWeights = reinterpret_cast<float *>(mxGetData(prhs[2]));
    //## Transformation Matrix
    float *transformationMatrix = reinterpret_cast<float *>(mxGetData(prhs[3]));
    //## Parameters
    //### Allow Scaling
    bool useScaling = static_cast<bool>(mxGetScalar(prhs[4]));
    
    //# Execute c++ function                                        
    meshmonk::compute_rigid_transformation_mex(floatingFeatures, numFloatingElements,
                                                correspondingFeatures, inlierWeights,
                                                transformationMatrix,
                                                useScaling);
    
    //# Set Output
//     plhs[0] = mxCreateNumericMatrix(numFloatingElements, 6, mxSINGLE_CLASS, mxREAL); // output: double matrix
//     auto output = mxGetPr(plhs[0]);
//     //## Copy result form c++ function into the output
//     for (unsigned i = 0 ; i < numFloatingElements * 6 ; i++){
//         output[i] = floatingFeatures[i];
//     }
  
}
