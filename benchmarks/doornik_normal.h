#ifndef __doornik_zignor__
#define	__doornik_zignor__
#include "../MT19937.h"
#include <math.h>

#define IRanU() rand_long64()
#define DRanU()	uniform_double_PRN()

static double DRanNormalTail(double dMin, int iNegative)
{
	double x, y;
	do
	{
        x = log(DRanU()) / dMin;
		y = log(DRanU());
    } while (-2 * y < x * x);
	return iNegative ? x - dMin : dMin - x;
}

#ifdef VIZIGNOR
#define ZIGNOR_C 256
#define ZIGNOR_R 3.6541528853609986
#define ZIGNOR_V 0.00492867323397
#define INT64_SCALE pow(2, -63)
#define REVERT_INT64_SCALE pow(2, 63)
#define REVERT_INT64_SCALE_SQUARED pow(2, 126)
#define MAX_INDEX 0xFF
#else
#define ZIGNOR_C 128 					/* number of blocks */
#define ZIGNOR_R 3.442619855899  /* start of the right tail */
					/* (R * phi(R) + Pr(X>=R)) * sqrt(2\pi) */
#define ZIGNOR_V 9.91256303526217e-3
#define MAX_INDEX 0x7F
#define REVERT_INT64_SCALE_SQUARED 1
#endif

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

#ifdef VIZIGNOR
    for (i = 0; i < iC; i++) 
    {
        s_adZigX[i] = INT64_SCALE*s_adZigX[i];
        s_adZigR[i] = REVERT_INT64_SCALE*s_adZigR[i];
    }
#endif
}

double DRanNormalZig(void)
{
	unsigned int i;
	double x, u, f0, f1;
	for (;;)
	{
#ifdef VIZIGNOR
		MT_FLUSH();
/*		Can save bits here, by using the leading 8 bits of the double random uniform in i */
        i = Rand->l & MAX_INDEX; 
		/* first try the rectangular boxes */
		if ((Rand->sl & 0x7fffffffffffffff) < s_adZigR[i])
			return Rand++->sl * s_adZigX[i];
		/* bottom box: sample from the tail */
		if (i == 0) 
            return DRanNormalTail(ZIGNOR_R, Rand++->sl < 0);
            /* is this a sample from the wedges? */
		x = Rand++->sl * s_adZigX[i];
#else
        u = 2 * DRanU() - 1;
		i = IRanU() & MAX_INDEX;
		/* first try the rectangular boxes */
		if (fabs(u) < s_adZigR[i])
			return u * s_adZigX[i];
		/* bottom box: sample from the tail */
		if (i == 0) 
            return DRanNormalTail(ZIGNOR_R, u < 0);
            /* is this a sample from the wedges? */
		x = u * s_adZigX[i];
#endif
		f0 = exp(-0.5 * (s_adZigX[i]   * s_adZigX[i]   * REVERT_INT64_SCALE_SQUARED - x * x) );
		f1 = exp(-0.5 * (s_adZigX[i+1] * s_adZigX[i+1] * REVERT_INT64_SCALE_SQUARED - x * x) );
		if (f1 + DRanU() * (f0 - f1) < 1.0)
			return x;
	}
}
#endif
