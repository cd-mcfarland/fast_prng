cdef extern from "math.h":
    long double erfl(long double x)

def erf(long double x): 
    return erfl(x)

