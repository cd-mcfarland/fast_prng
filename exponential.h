#include <stdlib.h>
#include <math.h>

#include "MT19937.h"

#include "exponential_layers.h"

void exponential_setup(){
	mt_init();
}

static inline double exponential_overhang(uint8_t i) {
  double x, y;
  dw128_t W = wide_uniform();
  if (W.d[0] < W.d[1]) { 
    x = W.d[0]; 
    y = 1 - W.d[1]; 
  }
  else { 
    x = W.d[1]; 
    y = 1 - W.d[0]; 
  }
  if (x < y - EXPONENTIAL_E[i]) return EXPONENTIAL_X[i] + EXPONENTIAL_dX[i]*x;
  x = EXPONENTIAL_X[i] + EXPONENTIAL_dX[i]*x;
  return EXPONENTIAL_Y[i] + EXPONENTIAL_dY[i]*y < exp(-x) ? x : exponential_overhang(i);
}

static inline double exponential() {
	MT_FLUSH();
	uint8_t i = Rand->s[0];
	if (i < EXPONENTIAL_BINS) {
	  Rand->l = (Rand->l >> 2) | EXP_SET;
	  return EXPONENTIAL_X[i]*(Rand++->d - 1);
	}
	i = Rand->i[1] > EXPONENTIAL_ipmf[Rand->s[1]] ? EXPONENTIAL_A[Rand->s[1]] : Rand->s[1]; 
	/* The last 8 bits of 64-bit random integer are beyond the resolution of a double precision float */
	Rand++;
	return i > 0 ? exponential_overhang(i) : EXPONENTIAL_X[0] + exponential();
}

