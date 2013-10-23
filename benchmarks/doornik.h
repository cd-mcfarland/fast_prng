#ifndef __doornik_zignor__
#define	__doornik_zignor__
#include "MT19937.h"
#include <math.h>

#define DRanU()	uniform_double()

static double DRanNormalTail(double dMin, int iNegative)
{
	double x, y;
	do
	{	x = log(DRanU()) / dMin;
		y = log(DRanU());
	} while (-2 * y < x * x);
	return iNegative ? x - dMin : dMin - x;
}

#define ZIGNOR_C 128 					/* number of blocks */
#define ZIGNOR_R 3.442619855899  /* start of the right tail */
					/* (R * phi(R) + Pr(X>=R)) * sqrt(2\pi) */
#define ZIGNOR_V 9.91256303526217e-3
 /* s_adZigX holds coordinates, such that each rectangle has*/
 /* same area; s_adZigR holds s_adZigX[i + 1] / s_adZigX[i] */
static double s_adZigX[ZIGNOR_C + 1], s_adZigR[ZIGNOR_C];

static void zigNorInit(int iC, double dR, double dV)
{
	int i; double f;

	f = exp(-0.5 * dR * dR);
	s_adZigX[0] = dV / f;  /* [0] is bottom block: V / f(R) */
	s_adZigX[1] = dR;
	s_adZigX[iC] = 0;

	for (i = 2; i < iC; ++i)
	{
		s_adZigX[i] = sqrt(-2 * log(dV / s_adZigX[i - 1] + f));
		f = exp(-0.5 * s_adZigX[i] * s_adZigX[i]);
	}
	for (i = 0; i < iC; ++i)
		s_adZigR[i] = s_adZigX[i + 1] / s_adZigX[i];
}

double DRanNormalZig(void)
{
	unsigned int i;
	double x, u, f0, f1;
	for (;;)
	{
/*		u = 2 * DRanU() - 1;
		i = IRanU() & 0x7F;
		Can save bits here, by using the leading 8 bits of the double random uniform in i
*/
		MT_FLUSH();
		i = Rand->s[0];
		Rand->l = (Rand->l & 0x000fffffffffffff) | 0x4000000000000000;
		u = Rand++->d - 3; 
		/* End of replacement code */

		/* first try the rectangular boxes */
		if (fabs(u) < s_adZigR[i])
			return u * s_adZigX[i];
		/* bottom box: sample from the tail */
		if (i == 0)
			return DRanNormalTail(ZIGNOR_R, u < 0);
		/* is this a sample from the wedges? */
		x = u * s_adZigX[i];
		f0 = exp(-0.5 * (s_adZigX[i] * s_adZigX[i] - x * x) );
		f1 = exp(-0.5 * (s_adZigX[i+1] * s_adZigX[i+1] - x * x) );
		if (f1 + DRanU() * (f0 - f1) < 1.0)
			return x;
	}
}
#endif
