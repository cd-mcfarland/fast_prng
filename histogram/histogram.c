#include <stdio.h>
#include <math.h>

#define TRIALS          	pow(10,9)
#define XMAX            	2
#define BINS            	2000
#define FILENAME		"counts.csv"

#ifdef EXPONENTIAL
#include "../exponential.h"
#define SETUP()			exponential_setup()
#define GENERATOR()		exponential()
#define NAME			"Exponential (Modified Ziggurat)"
#define XMIN            	0
#endif
#ifdef NORMAL
#include "../normal.h"
#define SETUP()			normal_setup()
#define GENERATOR()		normal()
#define NAME			"Standard Normal (Modified Ziggurat)"
#define XMIN            	-XMAX
#endif

int main( int argc, const char* argv[] ) {
    SETUP();

    // The intervals are the lower limit of each bin.
    double interval[BINS+1];
    int observed[BINS];
    for (int i = 0; i < BINS; i++) {
        interval[i] = XMIN + (XMAX - XMIN) * (i + 1.0) / BINS;
        observed[i] = 0;
    }
   
    // Create histogram
    for (int i = 0; i < TRIALS; i++) {
        double x = GENERATOR();
        if (x >= XMIN && x < XMAX) {
            	int index = (int)(BINS*(x - XMIN)/(XMAX - XMIN));
		observed[index]++;
	}
    }

    // Save the histogram
    FILE* f = fopen(FILENAME, "w");
    fprintf(f, "Low,High,%s\n", NAME);
    for (int i = 0; i < BINS; i++) {
	fprintf(f, "%.16e,%.16e,%d\n", interval[i], interval[i+1], observed[i]);
    }
    fclose(f);
}


