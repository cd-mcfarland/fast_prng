#ifndef __cdm_z_norm__
#define __cdm_z_norm__	
#include "MT19937.h"
#include "exponential.h"
#include "normal_layers.h"

#define N(x) (exp(-(x*x/2)))

void normal_setup() {
  mt_init();
}

static inline double norm_overhang(uint8_t i) {
  MT_FLUSH();
  dw128_t W;
  W.si = _mm_set_epi64x(Rand[0].l, Rand[1].l);	
  Rand	 += 2;
  W.si = _mm_or_si128(_mm_srli_epi64(W.si, 2), sse2_int_set); 
  W.sd = _mm_add_pd(W.sd, sse2_double_m_one);
  double x = norm_X[i] + (norm_X[i-1] - norm_X[i])*W.d[0], y = norm_Y[i] + (norm_Y[i-1] -norm_Y[i])*W.d[1];
  return y < N(x) ? x : norm_overhang(i);
}

static inline double norm_tail() {
  double x = exponential()/norm_X[0], y = exponential();
  return 2*y > x*x ? norm_X[0] + x : norm_tail();
}

static inline double normal() {
	MT_FLUSH();
	uint8_t i = Rand->s[0]; 
	if (i < NORM_BINS) {
      Rand->l = (Rand->l & 0x000fffffffffffff) | 0x4000000000000000;
	  return norm_X[i]*(Rand++->d - 3); 
	}	
/* Contingent upon Left to Right Associativity; I don't know if its good practice to write this in 1 line.*/

	i = Rand->i[1] >= norm_ipmf[Rand->s[1]] ? norm_A[Rand->s[1]] : Rand->s[1]; 
/* The last 8 bits of 64-bit random integer are beyond the resolution of a double precision float */
    if (i == 0) return (Rand++->s[2] > 127 ? 1 : -1)*norm_tail();
	return Rand++->s[2] > 127 ? norm_overhang(i) : -norm_overhang(i);
}
#endif

