from distutils.core import setup
from Cython.Build import cythonize

setup(
	name='cdm_ziggurat',
	version='1.0.0',
	description='Fast exponential and normal Pseudo Random Number Generator',
	author='Christopher McFarland',
	url='',
	maintainer_email='christopherdmcfarland+pypi@gmail.com',	
	classifiers =[
		'Programming Language :: Python',
		'Programming Language :: Python :: 3',
		'Programming Language :: Python :: 2.7',
		'Programming Language :: Cython',
		'Programming Language :: C',
		'License :: OSI Approved',
		'License :: OSI Approved :: MIT License',
		'Operating System :: Unix',
		'Intended Audience :: Science/Research',
		'Natural Language :: English',
		'Topic :: Scientific/Engineering',
		],
	license = 'MIT License',	
	packages = [''],
	ext_modules = cythonize("cdm_ziggurat.pyx")
)

