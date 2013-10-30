from numpy cimport *
import numpy as np

cdef extern from "../exponential.h":
	void exponential_setup()
	double c_exponential "exponential"()

exponential_setup()
def exponential(double scale=1.0, size=None):
	"""	exponential(scale=1.0, size=None)

Exponential distribution.

Its probability density function is

.. math:: f(x; \frac{1}{\beta}) = \frac{1}{\beta} \exp(-\frac{x}{\beta}),

for ``x > 0`` and 0 elsewhere. :math:`\beta` is the scale parameter,
which is the inverse of the rate parameter :math:`\lambda = 1/\beta`.
The rate parameter is an alternative, widely used parameterization
of the exponential distribution [3]_.

The exponential distribution is a continuous analogue of the
geometric distribution.  It describes many common situations, such as
the size of raindrops measured over many rainstorms [1]_, or the time
between page requests to Wikipedia [2]_.

Parameters
----------
scale : float
    The scale parameter, :math:`\beta = 1/\lambda`.
size : tuple of ints
    Number of samples to draw.  The output is shaped
    according to `size`.

See http://................................... for further details.
"""
	cdef long totalSize = np.multiply.reduce(size) if size else 1
	cdef double *element, *end	
	cdef ndarray[dtype=double, ndim=1] output = np.empty(totalSize)
	
	element = &(output[0]) 
	end = element + totalSize
		
	if scale != 1.0:
		while element < end:
			element[0] = scale*c_exponential()
			element += 1
	else:
		while element < end:
			element[0] = c_exponential()
			element += 1
	return output.reshape(size)


cdef extern from "../normal.h":
	void normal_setup()
	double c_normal "normal"()


normal_setup()
def normal(double loc=0.0, double scale=1.0, size=None):
	"""normal(loc=0.0, scale=1.0, size=None)

Draw random samples from a normal (Gaussian) distribution.

The probability density function:

		p(x) = 1/sqrt(2*pi*sigma**2)*exp( -(x - loc)**2/(2*sigma^2) ),

	where mu = loc, scale = sigma, and sqrt/pi/exp are defined in numpy. 

Parameters
----------
loc : float
    Mean ("centre") of the distribution.
scale : float
    Standard deviation (spread or "width") of the distribution.
size : tuple of ints
    Output shape.  If the given shape is, e.g., ``(m, n, k)``, then
    ``m * n * k`` samples are drawn.

see http://................................ for further details. 
"""
	cdef long totalSize = np.multiply.reduce(size) if size else 1
	cdef double *element, *end
	cdef ndarray[dtype=double, ndim=1] output = np.empty(totalSize)

	element = &(output[0])
	end = element + totalSize
	if scale != 1.0:
		if loc != 0:
			while element < end:
				element[0] = scale*c_normal() + loc
				element += 1
		else:
			while element < end:
				element[0] = scale*c_normal()
				element += 1
	else:
		if loc != 0:
			while element < end:
				element[0] = c_normal() + loc
				element += 1
		else:
			while element < end:
				element[0] = c_normal()
				element += 1
	return output.reshape(size)
