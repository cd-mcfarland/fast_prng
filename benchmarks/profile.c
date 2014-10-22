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
#define TRIALS 			pow(10, 9)
#include <stdio.h>

#ifdef FLOAT_EXPONENTIAL
#include "../float_exponential.h"
#define SETUP			exponential_setup()
#define GENERATOR()		exponential()
#define NAME			"32-bit exponential"
#endif
#ifdef EXPONENTIAL
#include "../exponential.h"
#define SETUP			exponential_setup()
#define GENERATOR()		exponential()
#define NAME			"exponential (Modified Ziggurat)"
#endif
#ifdef NORMAL
#include "../normal.h"
#define SETUP			normal_setup()
#define GENERATOR()		normal()
#define NAME			"Standard Normal (Modified Ziggurat)"
#endif
#ifdef POISSON
#include "../poisson.h"
#define SETUP			poisson_setup(1)
#define GENERATOR()		poisson()
#define NAME			"Poisson"
#endif
#ifdef DOORNIK
#include "doornik_normal.h"
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

int main(int argc, char *argv[]){
	double x = 0;
	double start_time = (double)clock();
	SETUP;
	int i;
	for (i=0; i<TRIALS; i++) x += (double)GENERATOR();
	double end_time = (double)clock();
	//Output
	printf("Created %ld %s distributed PRNs with mean %g in %g seconds.\n", (long)TRIALS, NAME, x/TRIALS, (end_time - start_time)/CLOCKS_PER_SEC);
}

