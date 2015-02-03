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


#include "./shared.h"
#include "../MT19937.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    mt_init();
    if (nrhs == 0) {
		plhs[0] = mxCreateDoubleScalar(uniform_double_PRN());
		return;		
	}
	
	plhs[0] = createOutputArray(nrhs, prhs);

	double *element = mxGetPr(plhs[0]);
	double *end = element + mxGetNumberOfElements(plhs[0]);
	
	while (element < end) {
		*element++ = uniform_double_PRN();
	}
}
