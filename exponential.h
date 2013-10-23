#ifndef __cdm_z_exp__
#define __cdm_z_exp__
#include <stdlib.h>
#include <math.h>

#include "MT19937.h"

#include "exponential_layers.h"
#define EXPONENTIAL_MINIMAL_TEST	0.092587156309399996

void exponential_setup(){
	mt_init();
}

static inline double exponential_overhang(uint8_t i) {
  double x, y;
  dw128_t W = wide_uniform();
  if (W.d[0] > W.d[1]) { 
    y = W.d[0];
    W.d[0] = W.d[1];
  }
  else {
    y = W.d[1];
  }
  x = EXPONENTIAL_X[i] + EXPONENTIAL_dX[i]*W.d[0];
  if ((y - W.d[0] <= EXPONENTIAL_MINIMAL_TEST) && (EXPONENTIAL_Y[i] + EXPONENTIAL_dY[i]*(1 - y) > exp(-x) )) return exponential_overhang(i);
  return x;
}

static inline double exponential() {
	MT_FLUSH();
	uint8_t i = Rand->s[0];
	if (i < EXPONENTIAL_BINS) {
	  Rand->l = (Rand->l >> 2) | EXP_SET;
	  return EXPONENTIAL_X[i]*(Rand++->d - 1);
	}
	Rand++;
	i = Rand->s[7];
	MT_FLUSH();
	if (Rand->i[0] > EXPONENTIAL_ipmf[i]) i = EXPONENTIAL_A[i];

//	i = Rand->i[1] >= EXPONENTIAL_ipmf[Rand->s[1]] ? EXPONENTIAL_A[Rand->s[1]] : Rand->s[1]; 
	// The last 8 bits of 64-bit random integer are beyond the resolution of a double precision float
	Rand++;
	return i == 0 ? EXPONENTIAL_X[0] + exponential() : exponential_overhang(i);
//	return i == 0 ? EXPONENTIAL_X[0] - log(uniform_double()) : exponential_overhang(i);
}
#endif

