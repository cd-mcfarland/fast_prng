#include <stdio.h>
#include <math.h>

#define TRIALS          pow(10,9)
#define XMAX            2
#define BINS            2000
#define FILENAME		"counts.csv"

#ifdef EXPONENTIAL
#include "../exponential.h"
#define SETUP()			exponential_setup()
#define GENERATOR()		exponential()
#define NAME			"Exponential (Modified Ziggurat)"
#define XMIN            0
#endif
#ifdef NORMAL
#include "../normal.h"
#define SETUP()			normal_setup()
#define GENERATOR()		normal()
#define NAME			"Standard Normal (Modified Ziggurat)"
#define XMIN            -XMAX
#endif

int main( int argc, const char* argv[] ) {
    SETUP();

        // Create histogram
    int observed[BINS] = {0};
    for (int i = 0; i < TRIALS; i++) {
        double x = GENERATOR();
        if (x >= XMIN && x < XMAX) {
            int index = (int)(BINS*(x - XMIN)/(XMAX - XMIN));
            observed[index]++;
        }

        // Save the histogram
    FILE* f = fopen(FILENAME, "w");
    fprintf(f, "Low,High,%s\n", NAME);
    for (int i = 0; i < BINS; i++) {
        double low = XMIN + (XMAX - XMIN) * i / BINS;
        double high = XMIN + (XMAX - XMIN) * (i+1) / BINS;
	    fprintf(f, "%.16e,%.16e,%d\n", low, high, observed[i]);
    }
    fclose(f);
}


