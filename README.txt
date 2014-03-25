OVERVIEW

    An unfinished manuscript can be found here:
    https://docs.google.com/document/d/1BRF8lhebODmYfku_tW0B38phyiOuplsU1UZuaOsRWfE/edit?usp=sharing


INSTALLATION & USAGE

    Usage in C/C++:

        Simply include this directory in the path of your compiler. See 
        benchmarks/profile.c for an example of usage. This program should 
        compile and profile the exponential algorithm with the command:

            gcc -O3 -DEXPONENTIAL profile.c -lm -o exponential_profiler.out

    Usage in Python:

        Install using Pip via https://pypi.python.org/pypi/fast_prng

    Usage in Matlab:

        Compile within the Matlab directory using: 

            mex cdm_exprnd.c
            mex cdm_randn.c

        Make sure the compiled functions are within the Matlab Path. Usage of 
        the two functions mimics exprnd and randn.      

CONTACT

    Feedback is most appreciated. If you have any questions or suggestions, 
    please contact me at christopherdmcfarland+fast_prng@gmail.com

WEBSITE

    All documentation and source code can be found at 
    https://bitbucket.org/cdmcfarland/fast_prng

FILE MANIFEST

    MT19937.h  

        Modified code of Super-Fast Mersene Twister used to generate uniform
        PRNs (http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/). Some code
        designed to support older architectures was removed, as this package
        was not tested on those architectures. Functions for seeding the PRNG 
        and accessing the array of PRNs were added.
    
    exponential.h  
    
        Exponential PRNG to be used with C/C++ code.         

    normal.h  

        Normal PRNG to be used with C/C++ code. 


    create_layers.py  

        Python script used to calculate X, i_max, A, f(X), and epsilon, as 
        described in the main text. All calculations are done to 128-bit 
        precision and then rounded to 64-bit precision to ensure no 
        rounding errors.
    
    erfl.pyx  

        Error function calculated to long double precision transfered to 
        python for create_layers.py
    
    quality_test.c 

        Returns raw moments of a set of PRNs. An example of usage is described
        in the Appendix of the main text.
    
    benchmarks/  
    
        Collection of scripts used to profile. 
        profile.c requires a definition at compile-time (EXPONENTIAL, NORMAL, 
        MARSAGLIA, DOORNIK) that chooses the algorithm to profile. E.g. usage:

          $ gcc -O3 -DNORMAL profile.c -lm -o normal.out
          $ ./normal.out 
          Created 1000000000 standard normal distributed pseudo-random numbers 
          with mean -3.18405e-05 in 3.33 seconds.

        Competing algorithms, modified to accept uniform PRNs from MT19937.h,
        are also provided. 
    
    matlab/  
    
default = https://cdmcfarland@bitbucket.org/cdmcfarland/scripts
        Extensions for Matlab/Octave.
