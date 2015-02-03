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

#include "./shared.h"
#include "../normal.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    normal_setup();
    if (nrhs == 0) {
        plhs[0] = mxCreateDoubleScalar(normal());
        return;
    }
	plhs[0] = createOutputArray(nrhs, prhs);
    
    double *element = mxGetPr(plhs[0]);
    double *end = element + mxGetNumberOfElements(plhs[0]);
    while (element < end) {
        *element++ = normal();
    }
}
