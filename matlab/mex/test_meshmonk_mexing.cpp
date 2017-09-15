#include "mex.h"
#include <meshmonk.hpp>

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
    if(nrhs != 5) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                      "5 inputs required.");
    }
    
    //# Get Inputs
    //float *targetFeatures = mxGetPr(prhs[1]); //matlab syntax
    float *floatingFeatures = reinterpret_cast<float *>(mxGetData(prhs[0]));
    mwSize nFloatingElements = mxGetM(prhs[0]);
    float *targetFeatures = reinterpret_cast<float *>(mxGetData(prhs[1])); //jiarui syntax
    mwSize nTargetElements = mxGetM(prhs[1]);
    float multiplier = static_cast<float>(mxGetScalar(prhs[4]));
    
    //# Execute c++ function
    meshmonk::test_meshmonk_mexing_raw(floatingFeatures, targetFeatures, nFloatingElements, nTargetElements, multiplier);
    
    //# Set Output
    plhs[0] = mxCreateDoubleMatrix(nFloatingElements, 6, mxREAL); // output: double matrix
    auto output = mxGetPr(plhs[0]);
    //## Copy result form c++ function into the output
    for (unsigned i = 0 ; i < nFloatingElements * 6 ; i++){
        output[i] = floatingFeatures[i];
    }
}