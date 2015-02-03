/* 
 * R = exprnd(mu)
 * R = exprnd(mu,m,n,...)
 * R = exprnd(mu,[m,n,...])
 *
 * Generates exponential random numbers with same behavior as native exprnd.
 *
 *
 * R = exprnd(mu) generates random numbers from the exponential distribution with mean parameter mu. 
 * mu can be a vector, a matrix, or a multidimensional array. The size of R is the size of mu.
 * 
 * R = exprnd(mu,m,n,...) or R = exprnd(mu,[m,n,...]) generates an m-by-n-by-... array containing 
 * random numbers from the exponential distribution with mean parameter mu. mu can be a scalar or an 
 * array of the same size as R.
 *
 * */

#include "./shared.h"
#include "../exponential.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	exponential_setup();
	if (nrhs == 0) {
		plhs[0] = mxCreateDoubleScalar(exponential());
		return;		
	}
	if (mxGetNumberOfElements(prhs[0]) != 1)
        mexErrMsgIdAndTxt("CDM_ZIGGURAT:exprnd:argtype", "mu must be scalar.");
    
    double mu = mxGetScalar(prhs[0]);
	if (nrhs == 1) {
		plhs[0] = mxCreateDoubleScalar(mu*exponential());
		return;
	}

	plhs[0] = createOutputArray(nrhs - 1, prhs + 1);

	double *element = mxGetPr(plhs[0]);
	double *end = element + mxGetNumberOfElements(plhs[0]);
	
	if (mu == 1) {
		while (element < end) 
            *element++ = exponential();
	}
	else { 
		while (element < end) 
            *element++ = mu*exponential();
	}
}
