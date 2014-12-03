/*
 *  randn Normally distributed pseudorandom numbers.
 *      R = randn(N) returns an N-by-N matrix containing pseudorandom values drawn
 *      from the standard normal distribution.  randn(M,N) or randn([M,N]) returns
 *      an M-by-N matrix. randn(M,N,P,...) or randn([M,N,P,...]) returns an
 *      M-by-N-by-P-by-... array. randn returns a scalar.  randn(SIZE(A)) returns
 *      an array the same size as A.
 *                       
 *      Note: The size inputs M, N, P, ... should be nonnegative integers.
 *      Negative integers are treated as 0.
 *                                
 *      The sequence of numbers produced by randn is determined by the settings of
 *      the uniform random number generator that underlies RAND, randn, and RANDI.
 *      randn uses one or more uniform random values to create each normal random
 *      value.  Control that shared random number generator using RNG.
 *                                                           
 *      Examples:
 *                                                                
 *      Example 1: Generate values from a normal distribution with mean 1
 *       and standard deviation 2.
 *         r = 1 + 2.*randn(100,1);
 *                                                                                          
 *      Example 2: Generate values from a bivariate normal distribution with
 *      specified mean vector and covariance matrix.
 *         mu = [1 2];
 *         Sigma = [1 .5; .5 2]; R = chol(Sigma);
 *         z = repmat(mu,100,1) + randn(100,2)*R; 	*/

#include <inttypes.h>
#include "mex.h"
#include "matrix.h"
#include "../normal.h"

#define ASSERT(a, b, c) { if (!(a)) mexErrMsgIdAndTxt(b, c); }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	ASSERT(nlhs == 1, "cdm_randn:maxlhs", "One output required.");
    normal_setup();
    if (nrhs == 0) {
        plhs[0] = mxCreateDoubleScalar(normal());
        return;
    }

    int64_t i, n_dims = nrhs == 1 ? mxGetNumberOfElements(prhs[0]) : nrhs; 
    mwSize dims[n_dims];

    if (nrhs == 1) {
        double *dims_p = mxGetPr(prhs[0]);
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
        *element++ = normal();
    }
}

