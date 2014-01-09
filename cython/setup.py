from distutils.core import setup, Extension
import os

Name = 'cdm_ziggurat'

# The easiest way to re-cythonize code: delete cdm_ziggurat.c (it is automatically generated by Cython)
USE_CYTHON = False if Name+'.c' in os.listdir() else True

ext = '.pyx' if USE_CYTHON else '.c'

extensions = [Extension(Name, [Name+ext], extra_compile_args=['-O3'] )]

if USE_CYTHON:
	from Cython.Build import cythonize
	extensions = cythonize(extensions)

setup(
	name=Name,
	version='1.0.0',
	description='Fast exponential and normal Pseudo Random Number Generator',
	author='Christopher McFarland',
	author_email='christopherdmcfarland+pypi@gmail.com',	
	maintainer_email='christopherdmcfarland+pypi@gmail.com',	
	url='',
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
	ext_modules = extensions 
)

