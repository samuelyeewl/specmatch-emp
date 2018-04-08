.. _installation:

Installation
============

Download code from

https://github.com/samuelyeewl/specmatch-emp

or clone git repo 

::

    $ git clone https://github.com/samuelyeewl/specmatch-emp

Then, simply run 

::

   $ python setup.py install

This downloads all dependencies and data files, and creates the
command line scripts. 

The python package setuptools is required to install SpecMatch-Emp.
If not already installed on your system, run

::

	$ pip install setuptools


Minimal Install
---------------

As of v0.3, it is possible to install SpecMatch-Emp without downloading
the full library. This will allow only the shifting code to be used.

::

	$ python setup.py install --shift-only


Software dependencies
---------------------

To install

- setuptools
- wget (facilitates downloading files).

Python Packages

To use SpecMatch-Emp:

- numpy
- scipy
- matplotlib
- h5py
- astropy (>=1.3)
- pandas (>=0.18)
- lmfit (>=0.9)

To build library:
- astroquery
- isochrones


Data Files
----------

These data files are automatically downloaded upon installation.

- ``library.h5``: the HDF5 archive of spectra.
- ``hires_telluric_mask.csv``: mask of telluric lines in HIRES spectra
- ``detrend.csv``: parameter calibration info
- ``nso.fits``, ``j26.532_adj.fits``, ``j72.718_adj.fits``,
  ``j59.1926_adj.fits`` spectra used as wavelength standards
 
