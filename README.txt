OVERVIEW

    This code is the implementation of a exponentially- and normally-distributed
    PseudoRandom Number Generator described at http://arxiv.org/abs/1403.6870

INSTALLATION & USAGE

    C/C++:

        Simply include this directory in the path of your compiler. For an 
        example of using the exponential PRNG library (exponential.h), see 
        benchmarks/profile.c. This program is compiled with the command:

            gcc -O2 -DEXPONENTIAL profile.c -lm -o exponential_profiler.out

        normal.h contains code for the normal PRNG.

    Python:

        Install using the Python Package Index:
            $ pip install fast_prng

    Usage in Matlab:

        Compile within the Matlab directory: 

            $ mex cdm_exprnd.c
            $ mex cdm_randn.c
            $ mex cdm_rand.c

        Make sure the compiled functions are within the Matlab Path. Usage of 
        the two functions mimics exprnd and randn.      

CONTACT

    Feedback is most appreciated. If you have any questions or suggestions, 
    please contact me at christopherdmcfarland [at] gmail [dot] com

WEBSITE

    All documentation and source code can be found at 
    https://bitbucket.org/cdmcfarland/fast_prng

FILE MANIFEST

    MT19937.h  

        Modified code of Super-Fast Mersene Twister used for uniform PRNG
        (http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/). Some code that
        supported older architectures was removed. Functions for seeding the 
        PRNG and accessing the array of PRNs without unnecessary bounds-
        checking were added.
    
    exponential.h  
    
        Exponential PRNG to be used with C/C++ code.         

    normal.h  

        Normal PRNG to be used with C/C++ code. 

    create_layers.py  

        Python script that calculates X, L_max, A, f(X), and epsilon (described
        in the main text). All calculations are done to 128-bit precision, and         
        then rounded to 64-bit precision to avoid rounding errors.
    
    erfl.pyx  

        128-bit precision Error function for create_layers.py
    
    quality_test.c 

        Returns raw moments of a set of PRNs. An example of usage is described
        in the Appendix of the main text.
    
    benchmarks/  
    
        Collection of scripts used to profile. 
        profile.c requires a definition at compile-time (EXPONENTIAL, NORMAL, 
        MARSAGLIA, DOORNIK) that chooses the algorithm to profile. E.g. usage:

          $ gcc -O2 -DNORMAL profile.c -lm -o normal.out
          $ ./normal.out
          Created 1000000000 Doornik Standard Normal distributed PRNs with mean 2.39616e-05.
          Startup time: 55 (us).
          Mean execution time (per PRN): 7.987 (ns).

        Competing algorithms are also provided. 
    
    matlab/  
    
        Extensions for Matlab/Octave.
