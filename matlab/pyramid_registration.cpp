#include "mex.h"
#include "/home/jonatan/projects/meshmonk/meshmonk/meshmonk.hpp"

/*
Note: if matlab complains about the std library, follow http://stackoverflow.com/a/41383926
and https://nl.mathworks.com/matlabcentral/answers/132527-in-mex-files-where-does-output-to-stdout-and-stderr-go
*/
class mystream : public std::streambuf //from https://nl.mathworks.com/matlabcentral/answers/132527-in-mex-files-where-does-output-to-stdout-and-stderr-go
{
protected:
virtual std::streamsize xsputn(const char *s, std::streamsize n) { mexPrintf("%.*s", n, s); return n; }
virtual int overflow(int c=EOF) { if (c != EOF) { mexPrintf("%.1s", &c); } return 1; }
};
class scoped_redirect_cout
{
public:
	scoped_redirect_cout() { old_buf = std::cout.rdbuf(); std::cout.rdbuf(&mout); }
	~scoped_redirect_cout() { std::cout.rdbuf(old_buf); }
private:
	mystream mout;
	std::streambuf *old_buf;
};

scoped_redirect_cout mycout_redirect;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    //# Check input
    //## Number of input arguments
    if(nlhs != 1) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs",
                      "One output required.");
    }
    //## Number of output arguments
    if(nrhs != 20) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                      "20 inputs required.");
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
    //### Number of pyramid layers
    mwSize numPyramidLayers = static_cast<mwSize>(mxGetScalar(prhs[7]));
    //### Starting downsample percentage for floating mesh
    float downsampleFloatStart = static_cast<float>(mxGetScalar(prhs[8]));
    //### Starting downsample percentage for target mesh
    float downsampleTargetStart = static_cast<float>(mxGetScalar(prhs[9]));
    //### Final downsample percentage for floating mesh
    float downsampleFloatEnd = static_cast<float>(mxGetScalar(prhs[10]));
    //### Final downsample percentage for target mesh
    float downsampleTargetEnd = static_cast<float>(mxGetScalar(prhs[11]));
    //### Use symmetric correspondences
    bool correspondencesSymmetric = static_cast<bool>(mxGetScalar(prhs[12]));
    //### Number of neighbours to use to compute corresponding points
    mwSize correspondencesNumNeighbours = static_cast<mwSize>(mxGetScalar(prhs[13]));
    //### Inlier kappa
    float inlierKappa = static_cast<float>(mxGetScalar(prhs[14]));
    //### Sigma of gaussian used in visco-elastic smoothing of deformation field
    float transformSigma = static_cast<float>(mxGetScalar(prhs[15]));
    //### Starting number of viscous smoothing iterations
    mwSize transformNumViscousIterationsStart = static_cast<mwSize>(mxGetScalar(prhs[16]));
    //### Final number of viscous smoothing iterations
    mwSize transformNumViscousIterationsEnd = static_cast<mwSize>(mxGetScalar(prhs[17]));
    //### Starting number of elastic smoothing iterations
    mwSize transformNumElasticIterationsStart = static_cast<mwSize>(mxGetScalar(prhs[18]));
    //### Final number of elastic smoothing iterations
    mwSize transformNumElasticIterationsEnd = static_cast<mwSize>(mxGetScalar(prhs[19]));
    
    //# Execute c++ function
//     meshmonk::test_meshmonk_mexing_raw(floatingFeatures, targetFeatures, numFloatingElements, numTargetElements, transformSigma);
    meshmonk::pyramid_registration_mex(floatingFeatures, targetFeatures,
                                numFloatingElements, numTargetElements,
                                floatingFaces, targetFaces,
                                numFloatingFaces, numTargetFaces,
                                floatingFlags, targetFlags,
                                numIterations, numPyramidLayers,
                                downsampleFloatStart, downsampleTargetStart,
                                downsampleFloatEnd, downsampleTargetEnd,
                                correspondencesSymmetric, correspondencesNumNeighbours,
                                inlierKappa,
                                transformSigma,
                                transformNumViscousIterationsStart, transformNumViscousIterationsEnd,
                                transformNumElasticIterationsStart, transformNumElasticIterationsEnd);
    
    //# Set Output
    plhs[0] = mxCreateDoubleMatrix(numFloatingElements, 6, mxREAL); // output: double matrix
    auto output = mxGetPr(plhs[0]);
    //## Copy result form c++ function into the output
    for (unsigned i = 0 ; i < numFloatingElements * 6 ; i++){
        output[i] = floatingFeatures[i];
    }
  
}