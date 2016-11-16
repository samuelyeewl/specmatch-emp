.. specmatch-emp documentation master file, created by
   sphinx-quickstart on Fri Aug  5 10:20:16 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to specmatch-emp's documentation!
=========================================

SpecMatch-Emp is a package to extract stellar parameters from their spectra, by matching the target spectrum against empirically observed spectra from stars with known properties.


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


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

