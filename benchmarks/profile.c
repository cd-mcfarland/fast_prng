/* 
 * Returns the raw moments of a pseudorandom number generator.
 * 
 * This program confirms that the distributions in this library behave as
 * intended. Unlike central moments, raw moments are always unbiased estimators 
 * of the expected value of the raw moment. (e.g. The sample variance is an 
 * unbiased estimator of variance only if adjusted by a factor of N/(N - 1).)
 *
 * */

#include <math.h>
#define TRIALS 			pow(10, 9)
#include <stdio.h>
#include "MT19937.h"

#ifdef EXPONENTIAL
#include "exponential.h"
#define SETUP			exponential_setup()
#define GENERATOR()		exponential()
#define NAME			"exponential"
#endif
#ifdef NORMAL
#include "normal.h"
#define SETUP			normal_setup()
#define GENERATOR()		normal()
#define NAME			"standard normal"
#endif
#ifdef POISSON
#include "poisson.h"
#define SETUP			poisson_setup(1)
#define GENERATOR()		poisson()
#define NAME			"Poisson"
#endif
#ifdef DOORNIK
#include "doornik.h"
#define NAME			"Doornik Standard Normal"
#define SETUP			mt_init(); zigNorInit(ZIGNOR_C, ZIGNOR_R, ZIGNOR_V);
#define GENERATOR()		DRanNormalZig()
#endif
#ifdef MARSAGLIA
#include "marsaglia_tsang_exponential.h"
#define NAME "Marsaglia & Tsang (2000) exponential"
static int ke[256];
static float fe[256], we[256];
#define SETUP			mt_init(); r4_exp_setup(ke, fe, we)
static unsigned long int *jsr_unused;
#define GENERATOR()		r4_exp(jsr_unused, ke, fe, we)
#endif

void main(int argc, char *argv[]){
	SETUP;
	double x = 0;
	int i;
	for (i=0; i<TRIALS; i++) x += (double)GENERATOR();

	//Output
	printf("Created %ld %s distributed pseudo-random numbers with mean %g \n", (long)TRIALS, NAME, x/TRIALS);
}
