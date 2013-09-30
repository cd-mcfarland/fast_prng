
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
#include "exponential.h"

#define ASSERT(a, b, c) { if (!(a)) mexErrMsgIdAndTxt(b, c); }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	ASSERT(nlhs == 1, "CDM_ZIGGURAT:exprnd:maxlhs", "Too many output arguments.");
	exponential_setup();
	if (nrhs >= 1) {
		ASSERT(mxGetNumberOfElements(prhs[0]) == 1, "CDM_ZIGGURAT:exprnd:argtype","mu must be scalar.");
        double mu = mxGetScalar(prhs[0]);
        if (nrhs >= 2) {
			if (nrhs == 2) {
				ASSERT(mxIsInt64(prhs[1]), "CDM_ZIGGURAT:exprnd:argtype", "Dimension array must be int64.");
				plhs[0] = mxCreateNumericArray(mxGetNumberOfElements(prhs[1]), mxGetData(prhs[1]), mxDOUBLE_CLASS, mxREAL);
			}
			else {
				int64_t i;
				mwSize *dims = mxMalloc(sizeof(mwSize)*nrhs-1);
				for (i=1; i<nrhs; i++) {
					ASSERT(mxIsInt64(prhs[i]), "CDM_ZIGGURAT:cxprnd:argtype", "Each dimension must be int64.");
					dims[i-1] = *mxGetData(prhs[i]);
				}
				plhs[0] = mxCreateNumericArray(nrhs-1, dims, mxDOUBLE_CLASS, mxREAL);
				mxFree(dims);
			}	
			double *element = mxGetPr(plhs[0]);
			double *end = element + mxGetNumberOfElements(plhs[0]);
			while (element < end) {
				*element = mu*exponential();
				element++;
			}
		}
		else {
			plhs[0] = mxCreateDoubleScalar(mu*exponential());
		}
    }
	else {
		plhs[0] = mxCreateDoubleScalar(exponential());
	}
}
