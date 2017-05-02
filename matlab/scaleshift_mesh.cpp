#include "mex.h"
#include "/home/jonatan/projects/meshmonk/meshmonk/meshmonk.hpp"
#include "mystream.cpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    //# Check input
    //## Number of input arguments
    if(nlhs != 1) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs",
                      "1 output required.");
    }
    //## Number of output arguments
    if(nrhs != 4) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                      "4 inputs required.");
    }
    
    //# Get Inputs
    //## Old Features
    float *oldFeatures = reinterpret_cast<float *>(mxGetData(prhs[0]));
    mwSize numOldElements = mxGetM(prhs[0]);
    std::cout << "Num Old Elements - " << numOldElements << std::endl;
    //## Old Original Indices
    int *oldOriginalIndices = reinterpret_cast<int *>(mxGetData(prhs[1]));
    //## New Features
    float *newFeatures = reinterpret_cast<float *>(mxGetData(prhs[2]));
    mwSize numNewElements = mxGetM(prhs[2]);
    std::cout << "Num New Elements - " << numNewElements << std::endl;
    //## New Original Indices
    int *newOriginalIndices = reinterpret_cast<int *>(mxGetData(prhs[3]));
    
    //# Execute c++ function
    meshmonk::scaleshift_mesh_mex(oldFeatures, numOldElements,
                                  oldOriginalIndices,
                                  newFeatures, numNewElements,
                                  newOriginalIndices);
    
    //# Set Output
    //## New Features
    plhs[0] = mxCreateNumericMatrix(numNewElements, 6, mxSINGLE_CLASS, mxREAL);
    auto output = mxGetPr(plhs[0]);
    //### Copy result from c++ function into the output
    for (unsigned i = 0 ; i < numNewElements * 6 ; i++){
        output[i] = newFeatures[i];
    }
  
}