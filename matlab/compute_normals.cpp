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
    if(nrhs != 3) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                      "3 inputs required.");
    }
    
    //# Get Inputs
    //## Positions
    float *positions = reinterpret_cast<float *>(mxGetData(prhs[0]));
    mwSize numElements = mxGetM(prhs[0]);
    //## Faces
    int *faces = reinterpret_cast<int *>(mxGetData(prhs[1]));
    mwSize numFaces = mxGetM(prhs[1]);
    //## Normals
    float *normals = reinterpret_cast<float *>(mxGetData(prhs[2]));
    
    
    //# Execute c++ function
    meshmonk::compute_normals_mex(positions, numElements,
                                faces, numFaces,
                                normals);
    
  
}