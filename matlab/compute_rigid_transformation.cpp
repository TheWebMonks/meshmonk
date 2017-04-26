#include "mex.h"
#include "/home/jonatan/projects/meshmonk/meshmonk/meshmonk.hpp"
#include "mystream.cpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    //# Check input
    //## Number of input arguments
    if(nlhs != 1) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs",
                      "One output required.");
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
    //## Corresponding Features
    float *correspondingFeatures = reinterpret_cast<float *>(mxGetData(prhs[1]));
    //## Inlier Weights
    float *inlierWeights = reinterpret_cast<float *>(mxGetData(prhs[2]));
    //## Parameters
    //### Allow Scaling
    bool allowScaling = static_cast<bool>(mxGetScalar(prhs[3]));
    std::cout << "Allow Scaling          - " << allowScaling << std::endl;
    
    //# Execute c++ function                                        
    meshmonk::compute_rigid_transformation_mex(floatingFeatures, numFloatingElements,
                                                correspondingFeatures, inlierWeights,
                                                allowScaling);
    
    //# Set Output
    plhs[0] = mxCreateDoubleMatrix(numFloatingElements, 6, mxREAL); // output: double matrix
    auto output = mxGetPr(plhs[0]);
    //## Copy result form c++ function into the output
    for (unsigned i = 0 ; i < numFloatingElements * 6 ; i++){
        output[i] = floatingFeatures[i];
    }
  
}