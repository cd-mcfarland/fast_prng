/* 
 * Returns the raw moments of a pseudorandom number generator.
 * 
 * This program confirms that the distributions in this library behave as
 * intended. Unlike central moments, raw moments are always unbiased estimators 
 * of the expected value of the raw moment. (e.g. The sample variance is an 
 * unbiased estimator of variance only if adjusted by a factor of N/(N - 1).)
 *
 * */

int factorial(int n) {
    if (n<=1) return(1);
    else n=n*factorial(n-1);
    return(n);
}

int double_factorial(int n) {
    if (n<=1) return(1);
    else n=n*double_factorial(n-2);
    return(n);
}

#define TRIALS 			pow(10, 10)
#define NUM_RAW_MOMENTS	6
#include <inttypes.h>
#include <stdio.h>

#ifdef EXPONENTIAL
#include "exponential.h"
#define SETUP			exponential_setup()
#define GENERATOR()		exponential()
#define NAME			"exponential"
#define EXPECTED(i)     (factorial(i))
#endif
#ifdef NORMAL
#include "normal.h"
#define SETUP			normal_setup()
#define GENERATOR()		normal()
#define NAME			"standard normal"
#define EXPECTED(i)     ( (i)%2 == 0 ? double_factorial((i)-1) : 0 )
#endif
#ifdef OLD_NORMAL
#include "old_normal.h"
#define SETUP			normal_setup()
#define GENERATOR()		normal()
#define NAME			"Old standard normal"
#define EXPECTED(i)     ( (i)%2 == 0 ? double_factorial((i)-1) : 0 )
#endif

int main(int argc, char *argv[]){
	SETUP;
	long i, j;
	double val, X[NUM_RAW_MOMENTS], x_j;
	for (i=0; i<NUM_RAW_MOMENTS; i++) {
		X[i] = 0;
	}

	for (i=0; i<TRIALS; i++) {
		val = (double)GENERATOR();
		for (j=0, x_j=val; j<NUM_RAW_MOMENTS; j++, x_j*=val) {
			X[j] += x_j;
		}
	}

	//Output
	printf("Created %ld %s distributed pseudo-random numbers...\n", (long)TRIALS, NAME);
	for (i=0; i<NUM_RAW_MOMENTS; i++) {
			printf("X%lu: %f (Expected %i)\n", i+1, X[i]/TRIALS, EXPECTED(i+1));	
	}
}
