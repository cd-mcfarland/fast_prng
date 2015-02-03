# Fast PRNG: an exponential- and normal PseudoRandom Number Generator

## OVERVIEW
-----------
This code is the implementation of a exponentially- and normally-distributed
PseudoRandom Number Generator described at [http://arxiv.org/abs/1403.6870].

## INSTALLATION & USAGE
-----------------------
### C/C++:

Source code and functions are found in exponential.h and normal.h. Simply 
include this directory in the path of your compiler. For an example of usage, 
see benchmarks/profile.c, which to test the exponential PRNG, can be compiled 
with the command:

    $ gcc -O2 -DEXPONENTIAL profile.c -lm -o exponential_profiler.out

### Python:

Use the Python Package Index:
            
    $ pip install fast_prng

The module's functions mimic numpy.random functions of the same name. 

### Matlab:

Ensure that the mex command is in the PATH, and then execute:
            
    $ matlab/install_matlab_functions.sh

This script will install to the default MATLAB userpath, unless $USERPATH is 
redefined. The installed functions: fast_exprnd, fast_randn, and fast_rand
behave identically to the native MATLAB functions that share the same name.  

##CONTACT
---------
Feedback is most appreciated. If you have any questions or suggestions, please 
contact me at christopherdmcfarland [at] gmail [dot] com

WEBSITE
-------
All documentation and source code can be found at 
[https://bitbucket.org/cdmcfarland/fast_prng]

## FILE MANIFEST
----------------
    
### MT19937.h  

Modified code of Super-Fast Mersene Twister used for uniform PRNG
[http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/]. Some code that 
supported older architectures was removed. Functions for seeding the PRNG and 
accessing PRNs without unnecessary bounds-checking were added.
    
### exponential.h  
    
Exponential PRNG to be used with C/C++ code.         

### normal.h  

Normal PRNG to be used with C/C++ code. 

### shared.h

Shared constants and MACROS used by exponential.h and normal.h.

### create_layers.py  

Python script that calculates X, L_max, A, f(X), and epsilon (described in the 
main text). All calculations are done to 128-bit precision, and then rounded to 
64-bit precision to avoid rounding errors.
    
### erfl.pyx  

128-bit precision Error function for create_layers.py
    
### quality_test.c 

Returns raw moments of a set of PRNs. An example of usage is described in the 
Appendix of the main text.
    
### benchmarks/  
    
A collection profiling scripts. profile.c requires a definition at compile-time 
(e.g. EXPONENTIAL, NORMAL) that selects the algorithm to profile. For example:

    $ gcc -O2 -DNORMAL profile.c -lm -o normal.out
    $ ./normal.out
    Created 1000000000 Doornik Standard Normal distributed PRNs with mean 2.39616e-05.
    Startup time: 55 (us).
    Mean execution time (per PRN): 7.987 (ns).

Competing algorithms are also provided. 
    
### matlab/  
    
Extensions for Matlab/Octave.
