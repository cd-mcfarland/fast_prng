
/*  CDM_RAND Uniformly distributed pseudorandom numbers.
 *
 *  R = CDM_RAND(N) returns an N-by-N matrix containing pseudorandom values i
 *  drawn from the standard uniform distribution on the open interval(0,1).  
 *  RAND(M,N) returns an M-by-N matrix.  RAND(M,N,P,...) or returns an 
 *  M-by-N-by-P-by-... array.  RAND returns a scalar.  
 *
 *  Note: The size inputs M, N, P, ... should be nonnegative integers.
 *  Negative integers are treated as 0.
 *                                
 *  The sequence of numbers produced by RAND is determined by the Super Fast
 *  Mersenne Twister Algorithm used for the other provided PRNGs. 
 */


#include <inttypes.h>
#include "mex.h"
#include "matrix.h"
#include "../MT19937.h"

#define ASSERT(a, b, c) { if (!(a)) mexErrMsgIdAndTxt(b, c); }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	ASSERT(nlhs == 1, "CDM_ZIGGURAT:rand:maxlhs", "One output required.");
    mt_init();
    if (nrhs == 0) {
		plhs[0] = mxCreateDoubleScalar(uniform_double_PRN());
		return;		
	}
	
	int64_t i, n_dims = nrhs == 1 ? mxGetNumberOfElements(prhs[1]) : nrhs;
	mwSize dims[n_dims];
	double *dims_p;

	if (nrhs == 1) {
		dims_p = mxGetPr(prhs[0]);
		for (i=0; i<n_dims; i++) {
			dims[i] = (int64_t)(*dims_p++);
		}
	}
	else {
		for (i=0; i<n_dims; i++) {
			dims[i] = (int64_t)(mxGetScalar(prhs[i]));
		}
	}

	plhs[0] = mxCreateNumericArray(n_dims, dims, mxDOUBLE_CLASS, mxREAL);

	double *element = mxGetPr(plhs[0]);
	double *end = element + mxGetNumberOfElements(plhs[0]);
	
	while (element < end) {
		*element = uniform_double_PRN();
		element++;
	}
}
