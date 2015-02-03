#pragma once
/* Shared includes and output array construction function. */

#include <inttypes.h>
#include "mex.h"
#include "matrix.h"

mxArray *createOutputArray(int n_dim_args, const mxArray *dim_args[]) {
	int64_t i, n_dims = n_dim_args == 1 ? mxGetNumberOfElements(dim_args[0]) : n_dim_args;
    mwSize dims[n_dims];
	if (n_dim_args == 1) {
        double *dims_p = mxGetPr(dim_args[0]);
		for (i=0; i<n_dims; i++) 
			dims[i] = (int64_t)(*dims_p++);
	}
	else {
		for (i=0; i<n_dims; i++)
			dims[i] = (int64_t)(mxGetScalar(*dim_args++));
	}
    return n_dims > 1 ? mxCreateNumericArray(n_dims, dims, mxDOUBLE_CLASS, mxREAL) : mxCreateDoubleMatrix(*dims, *dims, mxREAL);
} 
