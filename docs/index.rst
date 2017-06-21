.. specmatch-emp documentation master file, created by
   sphinx-quickstart on Fri Aug  5 10:20:16 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Empirical SpecMatch
===================

Empirical SpecMatch, ``SpecMatch-Emp`` for short, is an algorithm for
characterizing the properties of stars based on their optical
spectra. Target spectra are compared against a dense spectral library
of well-characterized touchstone stars. A key component of
``SpecMatch-Emp`` is the library of high resolution (R~55,000), high
signal-to-noise (> 100) spectra taken with Keck/HIRES by the
California Planet Search. Spectra in the library may be accessed by a
convenient Python API. For non-Python users spectra are also available
monolithic memory-efficient HDF5 archive and a la carte as individual
FITS files.

Code availability
~~~~~~~~~~~~~~~~~
The ``SpecMatch-Emp`` codebase is availble on GitHub. Simply clone the
repository and run ``python setup.py install``. For more details, go to
the :ref:`installation` page.

Basic Usage
~~~~~~~~~~~
The most basic method of using SpecMatch-Emp is to run the command line script:

::

   smemp specmatch rjXXX.XXXX.fits

This performs the shifting and matching process, prints the results and
saves it into a .h5 file which can be read by ``SpecMatch.read_hdf()`` in 
python.

For more details on using the command line interface, visit :ref:`cmdline`.
For a more general usage primer, check out the :ref:`quickstart` page.

Attribution
~~~~~~~~~~~

If you make use of the specmatch-emp spectral library or matching code,
please cite Yee et al. (2017) `[ADS Record] <http://adsabs.harvard.edu/abs/2017arXiv170100922Y>`_



Contents:

.. toctree::
   :maxdepth: 2

   installation
   quickstart
   cmdline
   build-library
   library
   spectrum
   specmatch
   match
   shift


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

