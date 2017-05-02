#include "mex.h"
#include "/home/jonatan/projects/meshmonk/meshmonk/meshmonk.hpp"
#include "mystream.cpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    //# Check input
    //## Number of input arguments
    if(nlhs != 0) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs",
                      "Zero LHS outputs required.");
    }
    //## Number of output arguments
    if(nrhs != 8) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                      "8 inputs required.");
    }
    
    //# Get Inputs
    //## Features
    float *features = reinterpret_cast<float *>(mxGetData(prhs[0]));
    mwSize numElements = mxGetM(prhs[0]);
    std::cout << "Num Elements - " << numElements << std::endl;
    //## Faces
    int *faces = reinterpret_cast<int *>(mxGetData(prhs[1]));
    mwSize numFaces = mxGetM(prhs[1]);
    std::cout << "Num Faces    - " << numFaces << std::endl;
    //## Flags
    float *flags = reinterpret_cast<float *>(mxGetData(prhs[2]));
    //## Downsampled Features (output)
    float *downsampledFeatures = reinterpret_cast<float *>(mxGetData(prhs[3]));
    mwSize numDownsampledElements = mxGetM(prhs[3]);
    std::cout << "Num Downsampled Elements - " << numDownsampledElements << std::endl;
    //## Downsampled Faces (output)
    int *downsampledFaces = reinterpret_cast<int *>(mxGetData(prhs[4]));
    mwSize numDownsampledFaces = mxGetM(prhs[4]);
    std::cout << "Num Downsampled Faces    - " << numDownsampledFaces << std::endl;
    //## Downsampled Flags (output)
    float *downsampledFlags = reinterpret_cast<float *>(mxGetData(prhs[5]));
    //## Original Indices
    int *originalIndices = reinterpret_cast<int *>(mxGetData(prhs[6]));

    //## Parameters
    //### Downsample Ratio (between 0.0 and 1.0)
    float downsampleRatio = static_cast<float>(mxGetScalar(prhs[7]));
    std::cout << "downsampleRatio                    - " << downsampleRatio << std::endl;
    
    //# Execute c++ function
    meshmonk::downsample_mesh_mex(features, numElements,
                                  faces, numFaces,
                                  flags,
                                  downsampledFeatures, numDownsampledElements,
                                  downsampledFaces, numDownsampledFaces,
                                  downsampledFlags,
                                  originalIndices,
                                  downsampleRatio);
    
//     //# Set Output
//     //## Downsampled Features
//     plhs[0] = mxCreateNumericMatrix(numDownsampledElements, 6, mxSINGLE_CLASS, mxREAL);
//     auto output0 = mxGetPr(plhs[0]);
//     //### Copy result from c++ function into the output
//     for (unsigned i = 0 ; i < numDownsampledElements * 6 ; i++){
//         output0[i] = downsampledFeatures[i];
//     }
//     
//     //## Downsampled Faces
//     plhs[1] = mxCreateNumericMatrix(numDownsampledElements, 3, mxUINT32_CLASS, mxREAL);
//     auto output1 = mxGetPr(plhs[1]);
//     //### Copy result from c++ function into the output
//     for (unsigned i = 0 ; i < numDownsampledFaces * 3 ; i++){
//         output1[i] = downsampledFaces[i];
//     }
//     
//     //## Downsampled Flags
//     plhs[2] = mxCreateNumericMatrix(numDownsampledElements, 1, mxSINGLE_CLASS, mxREAL);
//     auto output2 = mxGetPr(plhs[2]);
//     //### Copy result from c++ function into the output
//     for (unsigned i = 0 ; i < numDownsampledElements * 1 ; i++){
//         output2[i] = downsampledFlags[i];
//     }
//     
//     //## Original Indices
//     plhs[3] = mxCreateNumericMatrix(numDownsampledElements, 1, mxUINT32_CLASS, mxREAL);
//     auto output3 = mxGetPr(plhs[3]);
//     //### Copy result from c++ function into the output
//     for (unsigned i = 0 ; i < numDownsampledElements * 1 ; i++){
//         output3[i] = originalIndices[i];
//     }
    
  
}