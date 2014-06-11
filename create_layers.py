#!/usr/bin/python3
#
# This script creates the 3 pre-computed lookup tables: X, Y = f(X), and A, for both the exponential and normal PRNG.
# Values are calculated to longdouble precision and then rounded to double precision.
# Tables are output into a C header file.

import numpy as np
from numpy import abs

eps = 10*np.finfo(np.double).eps
size = 256

#TYPE = 'EXPONENTIAL' 
TYPE = 'NORMAL'

def check_equal(x,y): 
	"""check_equal(x, y) -> raise Assertion Error if values in x & y differ by more than double precision"""
	assert (np.abs(x - y,dtype=np.longdouble) < eps).all(), "{:}  {:}".format(x,y)

def fsolve(f, a, b, xtol=1e-18):
	""" fsolve(f, fprime, x0, xtol=1e-18) -> x: f(x) = 0; implemented to longdouble precision
		
	Uses Secant Method, requires a <= x < b there are no safeguards with this implementation.
"""
	assert a < b, "a must be lower bound, b must be higher bound"
	x = [b, a]
	F = list(map(f, x))
	for n in range(2, 100):
		x.append(x[n-1] - F[n-1]*(x[n-1] - x[n-2])/(F[n-1] - F[n-2]))
		if abs(x[n] - x[n-1]) < xtol: 
			return x[n]
		F.append(f(x[n]))
	return 0

def realign(pmf):
	"""realign(pmf) -> X, A: X & A allow random sampling from pmf in O(1) time.
		
		pmf: defines any arbitrary probability mass function on the domain [0, len(pmf)); needn't be normalized.
	
	A discrete random variable r, defined by pmf, is drawn from X, A as follows:
		r = numpy.random.randint(len(pmf))
		if X[r] > numpy.random.rand():
			r = A[r]
	
	See	http://scorevoting.net/WarrenSmithPages/homepage/sampling.abs for a description of this method. This is an implementation of the method to longdouble precision.
"""
	pmf = np.copy(pmf)
	pmf /= pmf.mean()
	L = len(pmf)
	A = np.arange(L)
	B = np.arange(L+1)
	X = np.r_[pmf, np.longdouble(2)]		 # X[L] = sentinel.
	i, j = 0, L
	while True:
		while X[B[i]]< 1:				# In ascending order, find x_(b_i) > 1.
			i += 1
		while X[B[j]] >= 1:			 	# In descending order, find x_(b_j) < 1.
			j -= 1
		if i >= j:						# If ascent passes descent, end
			break
		B[i], B[j] = B[j], B[i]			# Swap b_i, b_j
	i = j
	j += 1
										# At this point, X[B][:j] is < 1 and X[B][j:] is > 1 
										# Also, X[B[i]] is < 1. 
	while i>=0:
		while X[B[j]] <= 1:				# Find x_(b_j) that needs probability mass
			j += 1
		if j > L-1:						# Nobody needs probability mass, Done
			break						   
		X[B[j]] -= 1 - X[B[i]]			# Send all of x_(b_i) excess probability mass to x_(b_j)				(x_(b_i) is now done).
		A[B[i]] = B[j]
		if X[B[j]] < 1:					# If x_(b_j) now has too much probability mass, 
			B[i], B[j] = B[j], B[i]		# Swap, b_i, b_j, it becomes a donor.
			j += 1
		else:							# Otherwise, leave it as an acceptor
			i -= 1

	new_pmf = np.copy(X[:-1])
	for a_i, pmf_i in zip(A, X[:-1]):
		new_pmf[a_i] += 1 - pmf_i	
	check_equal(new_pmf, pmf)
	return X[:-1], A

# redefine Transcendental functions for proper precision
exp = lambda x: np.exp(x, dtype=np.longdouble)
sqrt = lambda x: np.sqrt(x, dtype=np.longdouble)
power = lambda x, y: np.power(x, y, dtype=np.longdouble)
from erfl import erf

oneHalf = 1/np.longdouble(2)

if TYPE == 'EXPONENTIAL':
	volume = 1/np.longdouble(size)
	f = lambda x: exp(-x)
	CDF = lambda x: 1 - f(x)
elif TYPE == 'NORMAL':
	f = lambda x: exp(-x*x*oneHalf)				# This function is more efficient to calculate than the normalized Gaussian Distribution
	volume = sqrt(2*np.pi)/np.longdouble(2*size)
	CDF = np.vectorize( lambda x: sqrt(np.pi/2)*erf(x/sqrt(2)) )

X = [fsolve(lambda x: x*f(x) - volume, 1, 10 if TYPE == 'EXPONENTIAL' else 4)]
Y = [f(X[0])]

while X[-1] != 0:
	X.append(fsolve(lambda x: x*(f(x) - Y[-1]) - volume, X[-1]/2, X[-1]))
	Y.append(f(X[-1]))

X = np.array(X)
dX = -np.diff(X)							# dx_i = x_i-1 - x_i; x_-1 = x: f(x) = 0 = inf
Y = np.array(Y)
dY = np.diff(Y)
check_equal(X[1:-1]*dY[:-1], volume)

V = -np.diff(CDF(np.r_[np.inf, X]))
V[1:] -= Y[:-1]*dX
V = np.r_[V, np.zeros(size - len(X-1))]

bins = len(X) - 1

check_equal(size - bins, V.sum()/volume) 

V /= V.mean()

pmf, Map = realign(V)

max_uint64 = power(2, 64)
max_int64 = power(2, 64 - 1)
max_uint56 = power(2, 56)

ipmf = np.uint64(pmf*np.longdouble(max_uint56))    # WDS sampler uses first 8 bits of 64-bit random uint to sample an 8-bit integer
ipmf[pmf >= 1] = max_uint56


if TYPE == 'EXPONENTIAL':
	m = dY/dX
	assert (m > Y[:-1]).all(),	'tangent line must be steeper than initial derivative'
	assert (m < Y[1: ]).all(),  'tangent line must be shallower than final derivative'
	E = (Y[1:]-m*(1-X[1:]-np.log(m)))/dY
	print('__EXP_MINIMAL_TEST__', np.uint64(E.max()*max_uint64))
	X /= max_int64 
	Y /= max_int64
elif TYPE == 'NORMAL':
	print('__NORM_TAIL_BEGIN__', X[0]) 
	X /= max_int64 
	Y /= max_int64

######### OUTPUT
output = open(TYPE.lower() + "_layers.h", 'w')

short_type = 'exp' if TYPE == 'EXPONENTIAL' else 'norm' 

BINS = "__{:}_BINS__".format(short_type.upper())
SIZE = "__{:}_SIZE__".format(short_type.upper())


output.write("""#define\t{BINS}\t{bins}
#define\t{SIZE}\t{size}
""".format(**locals()))


def writeDoubleArray(Str):
	assert len(eval(Str)) == bins + 1, 'improper array for this output' 
	output.write('static double __{short_type}_{Str}__[{BINS}+1] = {{ '.format(Str=Str, **globals()) + ', '.join(map(str, eval(Str))) + '};\n\n' )

for double_array in ['X', 'Y']:
	writeDoubleArray(double_array)

output.write(
"""static uint8_t __{short_type}_map__[{SIZE}] = {{ {:} }};

static {:}_t __{short_type}_ipmf__[{SIZE}] = {{ {:}u }};
""".format(', '.join(map(str, Map)), str(ipmf.dtype), 'u, '.join(map(str, ipmf)), **locals()))
output.close()

