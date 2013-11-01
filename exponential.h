#ifndef __cdm_z_exp__
#define __cdm_z_exp__
#include <stdlib.h>
#include <math.h>

#include "MT19937.h"

#include "exponential_layers.h"
#define EXPONENTIAL_MINIMAL_TEST	0.092587156309399996*pow(2, 64)
#define X_0							7.56927469415

void exponential_setup(){
	mt_init();
}

static inline double exponential_overhang(uint8_t i) {
	unsigned int y;
	double x;
	MT_FLUSH();

	if (Rand[0].l > Rand[1].l) { 
		y = Rand[0].l;
		Rand[0].l = Rand[1].l;
	}
	else {
		y = Rand[1].l;
	}
	x = EXPONENTIAL_X[i]*pow(2, 56) + (EXPONENTIAL_X[i-1] - EXPONENTIAL_X[i])*(Rand[0].l >> 8);
	Rand += 2;
	if ((y - Rand[-2].l <= EXPONENTIAL_MINIMAL_TEST) && (EXPONENTIAL_Y[i] + (EXPONENTIAL_Y[i-1] - EXPONENTIAL_Y[i])*(pow(2, 64) - y) > exp(-x) )) {
		return exponential_overhang(i);
	}
	return x;
}

static inline double exponential() {
	MT_FLUSH();
	uint8_t i = Rand->s[0];
	if (i < EXPONENTIAL_BINS) {
		return EXPONENTIAL_X[i]*(Rand++->l >> 8);
	}

    uint8_t j = Rand->i[0] >= EXPONENTIAL_ipmf[Rand->s[6]] ? EXPONENTIAL_A[Rand->s[6]] : Rand->s[6]; 
	/* The last 8 bits of 64-bit random integer are beyond the resolution of a double precision float */
	Rand++;
	return j == 0 ? X_0 + exponential() : exponential_overhang(j);
}
#endif

