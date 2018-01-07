#!/usr/bin/python3
import argparse
import numpy as np

## For longdouble Error Function ####
import pyximport; pyximport.install()
from erfl import erf
#####################################

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
description="""Creates the 4 pre-computed lookup tables for both the exponential and normal PRNG: 
    X = length of the ziggurat layers,
    Y = f(X), 
    map = probability table (for Alias Sampling),
    ipmf = Alias Table.
Values are calculated to longdouble precision and then rounded to double/single precision. Tables output to a C header file.""")

parser.add_argument('--type', default='exponential', choices={'exponential', 'normal'},help="PRNG type.")
parser.add_argument('--size', default=256, type=int, help="Number of Ziggurat layers.")
parser.add_argument('--precision', default='double', choices={'float', 'double'}, help='Precision of tables') 
args = parser.parse_args()

if args.precision == 'double':
    max_int = 0x7fffffffffffffff 
    max_uint = 0xffffffffffffffff 
    to_int = np.int64
    c_int_type = 'int64_t'
elif args.precision == 'float':
    max_int = 0x7fffffff
    max_uint = 0xffffffff
    to_int = np.int32
    c_int_type = 'int32_t'

# redefine Transcendental functions for proper precision
exp = lambda x: np.exp(x, dtype=np.longdouble)
sqrt = lambda x: np.sqrt(x, dtype=np.longdouble)
power = lambda x, y: np.power(x, y, dtype=np.longdouble)

oneHalf = 1/np.longdouble(2)

def check_equal(x,y): 
    """check_equal(x, y) -> raise Assertion Error if values in x & y differ by more than double precision"""
    assert (np.abs(x - y,dtype=np.longdouble) < np.finfo(np.double).eps).all(), "{:}  {:}".format(x,y)

def fsolve(f, a, b, xtol=1e-19, iterations=100, start_safe=False):
    """ fsolve(f, a, b xtol=1e-19) -> x: f(x) = 0; implemented to longdouble precision
        
    Tries Secant, then uses Bisection Method to find root.
    Requires a <= x < b.
"""
    assert a < b, "a must be lower bound, b must be higher bound"
    x = [b, a]
    F = list(map(f, x))
    if np.sign(F[0]) == np.sign(F[1]): 
        raise ValueError("there is no root within [a, b]")

    if not start_safe:
        for n in range(2, iterations):
            x.append(x[n-1] - F[n-1]*(x[n-1] - x[n-2])/(F[n-1] - F[n-2]))
            if np.abs(x[n] - x[n-1]) < xtol: 
                return x[n]
            F.append(f(x[n]))

    length = x[0] - x[1]
    iterations = int(np.ceil(np.log2(length/xtol)))
    for iteration in range(iterations): 
        c = x[1] + length
        F_c = f(c)
        if F_c == 0.:
            break
        i = int(np.sign(F_c) == np.sign(F[1]))
        x[i] = c
        F[i] = F_c
        length /= 2
    return c

def realign(pmf):
    """realign(pmf) -> X, A: X & A allow random sampling from pmf in O(1) time.
        
        pmf: defines any arbitrary probability mass function on the domain [0, len(pmf)); needn't be normalized.
    
    A discrete random variable r, defined by pmf, is drawn from X, A as follows:
        r = numpy.random.randint(len(pmf))
        if X[r] > numpy.random.rand():
            r = A[r]
    
    See    http://scorevoting.net/WarrenSmithPages/homepage/sampling.abs for a description of this method. This is an implementation of the method to longdouble precision.
"""
    pmf = np.copy(pmf)
    pmf /= pmf.mean()
    L = len(pmf)
    A = np.arange(L)
    B = np.arange(L+1)
    X = np.r_[pmf, np.longdouble(2)]         # X[L] = sentinel.
    i, j = 0, L
    while True:
        while X[B[i]]< 1:                # In ascending order, find x_(b_i) > 1.
            i += 1
        while X[B[j]] >= 1:                 # In descending order, find x_(b_j) < 1.
            j -= 1
        if i >= j:                        # If ascent passes descent, end
            break
        B[i], B[j] = B[j], B[i]            # Swap b_i, b_j
    i = j
    j += 1
                                        # At this point, X[B][:j] is < 1 and X[B][j:] is > 1 
                                        # Also, X[B[i]] is < 1. 
    while i>=0:
        while X[B[j]] <= 1:                # Find x_(b_j) that needs probability mass
            j += 1
        if j > L-1:                        # Nobody needs probability mass, Done
            break                           
        X[B[j]] -= 1 - X[B[i]]            # Send all of x_(b_i) excess probability mass to x_(b_j)                (x_(b_i) is now done).
        A[B[i]] = B[j]
        if X[B[j]] < 1:                    # If x_(b_j) now has too much probability mass, 
            B[i], B[j] = B[j], B[i]        # Swap, b_i, b_j, it becomes a donor.
            j += 1
        else:                            # Otherwise, leave it as an acceptor
            i -= 1

    new_pmf = np.copy(X[:-1])
    for a_i, pmf_i in zip(A, X[:-1]):
        new_pmf[a_i] += 1 - pmf_i    
    check_equal(new_pmf, pmf)
    return X[:-1], A

if args.type == 'exponential':
    volume = 1/np.longdouble(args.size)
    f = lambda x: exp(-x)
    CDF = lambda x: 1 - f(x)
elif args.type == 'normal':
    f = lambda x: exp(-oneHalf*x*x)                
    pi = np.longdouble(0)
    pi_str = '31415926535897932384626433832795'
    for i, digit in enumerate(pi_str): 
        pi += np.longdouble(digit)*power(10, -i)

    volume = sqrt(2*pi)/np.longdouble(2*args.size)
    alpha = sqrt(oneHalf*pi)
    CDF = np.vectorize( lambda x: alpha*erf(sqrt(oneHalf)*x) )

def generate_X_i(vol=volume, last_Y_i=0):
    lower_bound = 1
    upper_bound = 10 if args.type == 'exponential' else 4
    while lower_bound > vol:        # Obtaining another layer is clearly hopeless after this point
    # There are two solutions for X_i (a tall-skinny box and a long-flat box). We want the latter, 
    # so lower_bound is reduced gradually to avoid solving for the tall-skinny box. 
        try:
            X_i = fsolve(lambda x: x*(f(x) - last_Y_i) - vol, lower_bound, upper_bound)
            yield X_i
            last_Y_i = f(X_i)
            lower_bound = 0.9*X_i   
            upper_bound = X_i       
        except ValueError:
            lower_bound *= 0.9

X = np.array(list(generate_X_i(vol=volume)) + [0])
Y = np.vectorize(f)(X)

dX = -np.diff(X)                            # dx_i = x_i-1 - x_i; x_-1 = x: f(x) = 0 = inf
dY = np.diff(Y)
check_equal(X[1:-1]*dY[:-1], volume)

V = -np.diff(CDF(np.r_[np.inf, X]))
V[1:] -= Y[:-1]*dX

V_tail = V[0]/V.sum()
print("Fraction of remainder in tail: {:.4%}".format(V_tail))

V = np.r_[V, np.zeros(args.size - len(V))]

layers = len(X) - 1
check_equal(args.size - layers, V.sum()/volume) 

print("Mean V: {:g}".format((np.arange(len(V))*V).sum()/V.sum())) 
V /= V.mean()

pmf, Map = realign(V)

print("static {:} X_0 = {:};".format(args.precision, X[0]))

ipmf = to_int(pmf*max_uint - max_int)    # WDS sampler uses first 8 bits of random uint to sample an 8-bit integer
ipmf[pmf >= 1] = max_int 

__relevant_equalities__ = """
dY_i = Y_i+1 - Y_i
dX_i = X_i - X_i+1
-m_i = Y_i - Y_i+1 / X_i - X_i+1
e_i = -m_i-1(x - X_i) + F_i - e^-x 
0 =  -m_i-1 + e^-x
x_max = -log(m_i-1)
y_max = Y_i(x_max) - f(x_max)
de_max = y_max / dY
de_max = (+m_i-1*(+log(m_i-1) + X_i) + F_i - e^log(m_i-1) / dY
de_max = { m_i-1 * (log[m_i-1] + X_i) + F_i - m_i-1 } /dY
"""
m = dY/dX
if args.type == 'exponential':
    assert (m > Y[:-1]).all(), 'tangent line must be steeper than initial derivative'
    assert (m < Y[1: ]).all(), 'tangent line must be shallower than final derivative'
    epsilon = (Y[1:]-m*(1-X[1:]-np.log(m)))
    E = epsilon/dY
    print('static {:} iE_max = {:};'.format(c_int_type, to_int(E[3:].max()*max_int)))
    X /= max_int 
    Y /= max_int
elif args.type == 'normal':
    print('__NORM_TAIL_BEGIN__', X[0])
    inflection_point = 1
    i_inflection = (X > inflection_point).sum()
    assert X[i_inflection + 1] < inflection_point, "Inflection point lies exactly on boundary between two layers" 
    print('static uint_fast8_t i_inflection = ', i_inflection+1)

    E = np.zeros_like(X)
  
    for i in range(1, len(X)):
        delta_i = lambda x: x*f(x) - m[i-1]
        X_max = fsolve(delta_i, X[i], X[i-1] if i != i_inflection else inflection_point) # X_max[i_inflection] should be max_x f(x) - Y_i_inflection(x) (not max_x abs(f(x) - Y_i_inflection(x)) ), because we want to sample points in the concave region.
        assert X_max != 0, "Failed to find epsion_i for i = {:}".format(i)
        Y_i = lambda x: Y[i] - m[i-1]*(x - X[i])
        E[i] = np.abs( Y_i(X_max) - f(X_max) ) / dY[i-1] 

    iE = to_int(np.ceil(E*max_int))
    print('static {:} iE_max = {:}, iE_min = {:};'.format(c_int_type, iE[i_inflection:].max(), iE[:i_inflection].max()))

    X /= max_int 
    Y /= max_int

######### OUTPUT
output = open(args.type+"_layers.h", 'w')
short_type = 'exp' if args.type == 'exponential' else 'norm' 

output.write("\n#define  __{:}_LAYERS__  {:}\n".format(short_type.upper(), layers))

Arrays = dict(X=X, Y=Y, map=to_int(Map), ipmf=ipmf)
if args.type == 'normal':
    Arrays['iE'] = iE

for name, data in Arrays.items():
    output.write('static {dtype} __{short_type}_{name}__[{length}] = {{ {data} }};\n'.format(
        dtype=args.precision if data.dtype == np.float128 else (str(data.dtype) + '_t'),
        short_type=short_type,
        name=name,
        length=len(data),
        data=', '.join(map(str,data))))
    

output.close()
