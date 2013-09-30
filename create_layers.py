#!/usr/bin/python3.2
from numpy import *
from numpy import abs
eps = finfo(double).eps
size = 256

TYPE = 'EXPONENTIAL' 
#TYPE = 'NORMAL'

def check_equal(x,y): assert (abs(x - y,dtype=longdouble) < eps).all(), "{:}  {:}".format(x,y)

def fsolve(f, x0, fprime, xtol=1e-18):
	""" fsolve(f, fprime, x0, xtol=1e-18) -> x: f(x) = 0; implemented to longdouble precision
		
	Uses Newton's Method, then Secant Method if convergence fails. fprime and x0 (initial guess) are necessary; there are no safeguards with this implementation.
"""
	for n in range(100):	# first, try Newton's Method
		x1 = x0 - f(x0)/fprime(x0)
		if abs(x1-x0) < xtol: 
			return x1
		else: x0 = x1
	else:					# Newton's method failed to converge, now try Secant Method
		x = [x0,x0 - f(x0)/fprime(x0)]
		F = [f(x[0]), f(x[1])]
		for n in range(2,100):
			x.append(x[n-1] - F[n-1]*(x[n-1] - x[n-2])/(F[n-1] - F[n-2]))
			if abs(x[n]-x[n-1]) < xtol: 
				return x[n]
			F.append(f(x[n]))
		else: 
			return 0

def realign(pmf):
	"""realign(pmf) -> X, A: X & A allow random sampling from pmf in O(1) time.
		
		pmf: defines any arbitrary probability mass function on the domain [0, len(pmf)).
	
	A discrete random variable r, defined by pmf, is drawn from X, A as follows:
		r = numpy.random.randint(len(pmf))
		if X[r] > numpy.random.rand():
			r = A[r]
	
	See	http://scorevoting.net/WarrenSmithPages/homepage/sampling.abs for a description of this method. This is an implementation of the method to longdouble precision.
"""
	L = len(pmf)
	A, B, X = arange(L), arange(L+1), r_[pmf,longdouble(2)]		 # X[L] = sentinel.
	i, j = 0, L
	while True:
		while X[B[i]]< 1:					   #In ascending order, find x_(b_i) > 1.
			i += 1
		while X[B[j]] >= 1:			 #In descending order, find x_(b_j) < 1.
			j -= 1
		if i >= j:									  #If ascent passes descent, end
			break
		B[i], B[j] = B[j], B[i]		 # Swap b_i, b_j
	i = j
	j += 1
		# At this point, X[B][:j] is < 1 and X[B][j:] is > 1 
		# Also, X[B[i]] is < 1. 
	while i>=0:
		while X[B[j]] <= 1:			 # Find x_(b_j) that needs probability mass
			j += 1
		if j > L-1:									 # Nobody needs probability mass, Done
			break								   
		X[B[j]] -= 1 - X[B[i]]		  # Send all of x_(b_i) excess probability mass to x_(b_j)				(x_(b_i) is now done).
		A[B[i]] = B[j]
		if X[B[j]] < 1:				 # If x_(b_j) now has too much probability mass, 
			B[i], B[j] = B[j], B[i] # Swap, b_i, b_j, it becomes a donor.
			j += 1
		else:										   # Otherwise, leave it as an acceptor
			i -= 1
	return X[:-1], A

if TYPE == 'EXPONENTIAL':
	f = lambda x: size*exp(-x,dtype=longdouble)
elif TYPE == 'NORMAL':
	heigth = exp(1/longdouble(2), dtype=longdouble)/sqrt(2*pi,dtype=longdouble)
	f = lambda x: size*heigth*exp(-x*x,dtype=longdouble)

X = []
Y = [0]

for i in range(size):
	X.append(fsolve(lambda x_i: x_i*(f(x_i) - Y[i]) - 1, X[i-1] if i!=0 else 10, fprime=lambda x_i: f(x_i)*(1-x_i) - Y[i]))
	Y.append(f(X[i]))
	if X[i] == 0:
		break
	check_equal(X[i]*(Y[i+1] - Y[i]), 1)
X = array(X)
bins = len(X) - 1
Y = array(Y)

dX = r_[longdouble(999), X[:-1] - X[1:]]								# dx_i = x_i-1 - x_i; x_-1 = x: f(x) = 0 = inf
V = r_[Y[1:] - Y[:-1]*(1+dX), zeros(size - len(X), dtype=longdouble)]	# len(V) must equal size

check_equal(size, V.sum() + bins)
Y /= size										  	# normalize back to 1.
V *= size/V.sum()								  	# normalize to size.
dY = Y[1:] - Y[:-1]									# dy_i = y_i+1 - y_i

pmf, A = realign(V)
newV = copy(pmf)
for A_i, pmf_i in zip(A, pmf):
	newV[A_i] += 1-pmf_i
for V_i, newV_i in zip(V, newV):
	check_equal(V_i, newV_i) 


ipmf = uint32(pmf*longdouble(2**32))
ipmf[pmf >= 1] = -1

output = open(TYPE.lower() + "_constants.h", 'w')


output.write("""#define\t{TYPE}_BINS\t{bins}
#define\t{TYPE}_SIZE\t{size}
""".format(**locals()))

def writeDoubleArray(Str):
	assert len(eval(Str)) == bins + 1, 'improper array for this output' 
	output.write('static double {TYPE}_{Str}[{TYPE}_BINS+1] = {{ '.format(Str=Str, **globals()) + ', '.join(map(str, eval(Str))) + '};\n\n' )

if TYPE == 'EXPONENTIAL':
	m = -dY/dX
	for i, (ym_i, y_i, y_ip1) in enumerate(zip(-m, Y[:-1], Y[1:])):
		assert ym_i > y_i and ym_i < y_ip1, 'i: {:} ym_i: {:} y_i: {:} y_i+1: {:}'.format(i, ym_i, y_i, y_ip1)
	E = (Y[1:]+m*(1-X-log(-m,dtype=longdouble)))/dY
	E1 = 1 - E
	writeDoubleArray('E')
	writeDoubleArray('E1')

Y = Y[:-1]

for double_array in ['X', 'Y', 'dX', 'dY']:
	writeDoubleArray(double_array)

output.write(
"""static uint8_t {TYPE}_A[{TYPE}_SIZE] = {{ {:} }};

static {:}_t {TYPE}_ipmf[{TYPE}_SIZE] = {{ {:}u }};
""".format(', '.join(map(str, A)), str(ipmf.dtype), 'u, '.join(map(str, ipmf)), **locals()))
output.close()
for i, x_i, dx_i, y_i, dy_i in zip(range(len(X)), X, dX, Y, dY):
	if i==0: 
		continue
	check_equal(f(x_i)/size, y_i + dy_i)
	check_equal(y_i, f(x_i + dx_i)/256)

