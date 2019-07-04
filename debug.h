#ifndef __debug__
#define __debug__

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#define I(x) { printf("%s: %ld\n",#x, (long)x); }
#define F(x) { printf("%s: %g\n",#x,x); }

#define VASSERT(condition) { if(!(condition)) { printf("ASSERT FAILED: %s @ %s (line %d)\n", #condition, __FILE__,  __LINE__ ); exit(EXIT_FAILURE);} }

static double __debug_sum__ = 0;
static long __debug_N__ = 0;
static char *__debug_name__;
static double __debug_min__ = INFINITY;
static double __debug_max__ = -INFINITY;

void _count(double x, char *name) {
	__debug_name__ = name;
	__debug_sum__ += x;	
	__debug_N__++;
    __debug_min__ = __debug_min__ < x ? __debug_min__ : x;
    __debug_max__ = __debug_max__ > x ? __debug_max__ : x;
}

void _describe(void) {
	static int firstTime = 1; 
	if (firstTime) {
		printf("Mean of %s is %g\n", __debug_name__, __debug_sum__/(double)__debug_N__);
		printf("Min of %s is %g\n", __debug_name__, __debug_min__);
		printf("Max of %s is %g\n", __debug_name__, __debug_max__);
		firstTime = 0;
	}
}

#define DESCRIBE(x) {_count(x, #x); atexit(_describe); }

#endif
