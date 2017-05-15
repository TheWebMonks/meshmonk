#include "mex.h"
#include <meshmonk/meshmonk.hpp>
#include "mystream.cpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    //# Check input
    //## Number of input arguments
    if(nlhs != 0) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs",
                      "Zero left-hand-side outputs required.");
    }
    //## Number of output arguments
    if(nrhs != 9) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                      "9 right-hand-side inputs required.");
    }
    
    //# Get Inputs
    //## Floating Features
    float *floatingFeatures = reinterpret_cast<float *>(mxGetData(prhs[0]));
    mwSize numFloatingElements = mxGetM(prhs[0]);
    //## Target Features
    float *targetFeatures = reinterpret_cast<float *>(mxGetData(prhs[1]));
    mwSize numTargetElements = mxGetM(prhs[1]);
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
    //### Number of neighbours to use to compute corresponding points
    mwSize correspondencesNumNeighbours = static_cast<mwSize>(mxGetScalar(prhs[7]));
    //### Flag threshold to mark corresponding flag as 0.0 or 1.0
    float correspondencesFlagThreshold = static_cast<float>(mxGetScalar(prhs[8]));
    
    //# Execute c++ function
    meshmonk::compute_correspondences_mex(floatingFeatures, targetFeatures,
                                        numFloatingElements, numTargetElements,
                                        floatingFlags, targetFlags,
                                        correspondingFeatures, correspondingFlags,
                                        correspondencesSymmetric, correspondencesNumNeighbours,
                                        correspondencesFlagThreshold);
    
//     //# Set Output
//     //## Corresponding Features
//     mwSize numCols = 6;
//     plhs[0] = mxCreateNumericMatrix(numFloatingElements, numCols, mxSINGLE_CLASS, mxREAL); // output: double matrix
//     auto output0 = mxGetPr(plhs[0]);
//     //## Copy result form c++ function into the output
//     for (unsigned i = 0 ; i < numFloatingElements * numCols ; i++){
//         output0[i] = correspondingFeatures[i];
//     }
//     
//     std::cout << "Finished copying corresponding Features " << std::endl;
//     
//     //## Corresponding Flags
//     numCols = 1;
//     plhs[1] = mxCreateNumericMatrix(numFloatingElements, numCols, mxSINGLE_CLASS, mxREAL);
//     std::cout << "plhs[1] : " << plhs << std::endl;
//     auto output1 = mxGetPr(plhs[1]);
//     //## Copy result form c++ function into the output
//     for (unsigned i = 0 ; i < numFloatingElements * numCols ; i++){
//         output1[i] = correspondingFlags[i];
//         std::cout << i << " - " << correspondingFlags[i] << std::endl;
//     }
//     
//     std::cout << "Finished copying corresponding flags " << std::endl;
  
}