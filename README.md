# Accelerator System Simulation Engine

This code base includes:

* Back-end physics simulator (in C),
* Top-level code to configure and run simulations (in Python),
* Logic to glue C with Python (SWIG).

Full documenation can be found in the `doc/` directory.

## Dependencies

* Python 2.7
* Numpy, Scipy, Matplotlib - Python packages. The code was developed with the following versions:

		In [1]: import numpy, scipy, matplotlib
		In [2]: numpy.__version__
		Out[2]: '1.8.2'
		In [3]: scipy.__version__
		Out[3]: '0.13.3'
		In [4]: matplotlib.__version__
		Out[4]: '1.3.1'

* SWIG:

The version used for development:

		$ swig -version 
		SWIG Version 3.0.5
		Compiled with g++ [x86_64-unknown-linux-gnu]
		Configured options: +pcre

This version has also been tested to work:

		$ swig -version
		SWIG Version 2.0.7
		Compiled with g++ [x86_64-unknown-linux-gnu]
		Configured options: +pcre

* Pydot - Python wrapper of Graphviz, used for generating images of connectivities.
* [Oct2py](https://pypi.python.org/pypi/oct2py): Python to GNU Octave bridge --> run m-files from Python for some unit tests.

## Dependencies

To run a basic run of unit tests, type (from the top directory):
		
		$ make unit_tests.log