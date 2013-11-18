
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

#include <inttypes.h>
#include "mex.h"
#include "matrix.h"
#include "../exponential.h"

#define ASSERT(a, b, c) { if (!(a)) mexErrMsgIdAndTxt(b, c); }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	ASSERT(nlhs == 1, "CDM_ZIGGURAT:exprnd:maxlhs", "One output required.");
	exponential_setup();
	if (nrhs == 0) {
		plhs[0] = mxCreateDoubleScalar(exponential());
		return;		
	}
	ASSERT(mxGetNumberOfElements(prhs[0]) == 1, "CDM_ZIGGURAT:exprnd:argtype","mu must be scalar.");
        double mu = mxGetScalar(prhs[0]);
	if (nrhs == 1) {
		plhs[0] = mxCreateDoubleScalar(mu*exponential());
		return;
	}
	
	int64_t i, n_dims = nrhs == 2 ? mxGetNumberOfElements(prhs[1]) : nrhs - 1;
	mwSize dims[n_dims];
	double *dims_p;

	if (nrhs == 2) {
		dims_p = mxGetPr(prhs[1]);
		for (i=0; i<n_dims; i++) {
			dims[i] = (int64_t)(*dims_p++);
		}
	}
	else {
		for (i=0; i<n_dims; i++) {
			dims[i] = (int64_t)(mxGetScalar(prhs[i+1]));
		}
	}
	plhs[0] = mxCreateNumericArray(n_dims, dims, mxDOUBLE_CLASS, mxREAL);
	double *element = mxGetPr(plhs[0]);
	double *end = element + mxGetNumberOfElements(plhs[0]);
	while (element < end) {
		*element = mu*exponential();
		element++;
	}
}
