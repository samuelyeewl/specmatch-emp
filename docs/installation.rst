.. _installation:

Installation
============

Simply run 

::

   $ python setup.py install

This downloads all dependencies and data files, and creates the
command line scripts. 

Software dependencies
---------------------

To install

- wget (facilitates downloading files).

Python Packages

To use SpecMatch-Emp:

- numpy
- scipy
- matplotlib
- h5py
- pandas (>=0.18)
- lmfit (>=0.9)

To build library:

- astropy
- astroquery
- isochrones


Data Files
~~~~~~~~~~

These data files are automatically downloaded upon installation.

- ``library.h5``: the HDF5 archive of spectra.
- ``hires_telluric_mask.csv``: mask of telluric lines in HIRES spectra
- ``detrend.csv``: parameter calibration info
- ``nso.fits``, ``rj26.532.fits``, ``rj72.718.fits``,
  ``rj59.1926.fits`` spectra used as wavelength standards
 
