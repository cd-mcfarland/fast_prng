/* 
 * Returns the execution time of a PseudoRandom Number Generator (PRNG). A mean
 * is also provided as a sanity check. 
 * 
 * This program is used to profile the various Ziggurat Algorithms available 
 * today. It profiles algorithms provided in this package along with algorithms 
 * provided in other recent publications. You must define the algorithm to 
 * profile during compilation, e.g.
 *
 * gcc -O2 -DEXPONENTIAL -o profile.out profile.c 
 *
 * */
#include <math.h>
#define TRIALS          pow(10,9)
#include <stdio.h>

#ifdef MT19937
#include "../MT19937.h"
#define SETUP()         mt_init()
#define GENERATOR()     uniform_double_PRN()
#define NAME            "64-bit random unsigned ints (sf-MT)"
#endif
#ifdef FLOAT_EXPONENTIAL
#include "../float_exponential.h"
#define SETUP()			exponential_setup()
#define GENERATOR()		exponential()
#define NAME			"32-bit exponential"
#endif
#ifdef EXPONENTIAL
#include "../exponential.h"
#define SETUP()			exponential_setup()
#define GENERATOR()		exponential()
#define NAME			"exponential (Modified Ziggurat)"
#endif
#ifdef NORMAL
#include "../normal.h"
#define SETUP()			normal_setup()
#define GENERATOR()		normal()
#define NAME			"Standard Normal (Modified Ziggurat)"
#endif
#ifdef DOORNIK
#include "doornik_normal.h"
#define NAME			"Doornik Standard Normal"
#define SETUP()			mt_init(); zigNorInit(ZIGNOR_C, ZIGNOR_R, ZIGNOR_V);
#define GENERATOR()		DRanNormalZig()
#endif
#ifdef MARSAGLIA
#include "marsaglia_tsang_exponential.h"
#ifdef DROP_ASSIGNMENT
#define NAME "Marsaglia (without x assignment)"
#else
#define NAME "Marsaglia & Tsang (2000) exponential"
#endif
static int ke[256];
static float fe[256], we[256];
#define SETUP()			mt_init(); r4_exp_setup(ke, fe, we);
static unsigned long int *jsr_unused;
#define GENERATOR()		r4_exp(jsr_unused, ke, fe, we)
#endif
#ifdef DOUBLE_MARSAGLIA
#include "double_marsaglia.h"
#define NAME "Double precision Marsaglia (additional options)"
#define SETUP()			mt_init(); r8_exp_setup();
#define GENERATOR()		r8_exp()
#endif

int main(int argc, char *argv[]){
	double x = 0;
	double start_time = (double)clock();
	SETUP();
	long i;
    double execution_start_time = (double)clock();
#pragma omp parallel for reduction(+:x) private(i)
    for (i=0; i<TRIALS; ++i) x += (double)GENERATOR();
	double end_time = (double)clock();
	//Output
	printf("Created %ld %s distributed PRNs with mean %g.\n", (long)TRIALS, NAME, x/TRIALS);
    printf("Startup time: %g (us).\n", (execution_start_time - start_time)/CLOCKS_PER_SEC*1e6);
    printf("Mean execution time (per PRN): %g (ns).\n", (end_time - execution_start_time)/CLOCKS_PER_SEC/TRIALS*1e9);
    return 0;
}

