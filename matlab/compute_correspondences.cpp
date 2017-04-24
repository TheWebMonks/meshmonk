#include "mex.h"
#include "/home/jonatan/projects/meshmonk/meshmonk/meshmonk.hpp"
#include "mystream.cpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    //# Check input
    //## Number of input arguments
    if(nlhs != 2) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs",
                      "One output required.");
    }
    //## Number of output arguments
    if(nrhs != 8) {
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
    //## FLoating Flags
    float *floatingFlags = reinterpret_cast<float *>(mxGetData(prhs[2]));
    //## Target Flags
    float *targetFlags = reinterpret_cast<float *>(mxGetData(prhs[3]));
    //## Corresponding Features (output)
    float *correspondingFeatures = reinterpret_cast<float *>(mxGetData(prhs[4]));
    //## Corresponding Flags (output)
    float *correspondingFlags = reinterpret_cast<float *>(mxGetData(prhs[5]));
    //## Parameters
    //### Use symmetric correspondences
    bool correspondencesSymmetric = static_cast<bool>(mxGetScalar(prhs[6]));
    std::cout << "Symmetric Correspondences          - " << correspondencesSymmetric << std::endl;
    //### Number of neighbours to use to compute corresponding points
    mwSize correspondencesNumNeighbours = static_cast<mwSize>(mxGetScalar(prhs[7]));
    std::cout << "Num Neighbours                     - " << correspondencesNumNeighbours << std::endl;
    
    //# Execute c++ function
    meshmonk::compute_correspondences_mex(floatingFeatures, targetFeatures,
                                        numFloatingElements, numTargetElements,
                                        floatingFlags, targetFlags,
                                        correspondingFeatures, correspondingFlags,
                                        correspondencesSymmetric, correspondencesNumNeighbours);
    
    //# Set Output
    //## Corresponding Features
    plhs[0] = mxCreateDoubleMatrix(numFloatingElements, 6, mxREAL); // output: double matrix
    auto output0 = mxGetPr(plhs[0]);
    //## Copy result form c++ function into the output
    for (unsigned i = 0 ; i < numFloatingElements * 6 ; i++){
        output0[i] = correspondingFeatures[i];
    }
    //## Corresponding Flags
    plhs[1] = mxCreateDoubleMatrix(numFloatingElements, 1, mxREAL); // output: double matrix
    auto output1 = mxGetPr(plhs[1]);
    //## Copy result form c++ function into the output
    for (unsigned i = 0 ; i < numFloatingElements ; i++){
        output1[i] = correspondingFlags[i];
    }
  
}