Building the Library
====================

The SpecMatch libary is stored as a monolithic HDF5 file and is
downloaded upon installation of the code. It's is constructed from a
number of different data files. For simplicity, we do not make all the
source files public. This document serves as a cookbook for CPS.

Establish Wavelength References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Raw spectra from CPS are not wavelength-calibrated. They are stored as
a list of intensities and pixel values. Spectra are resampled onto the
rest-wavelength scale using the NSO solar spectrum by
cross-correlating segments of spectrum. Cross-correlating spectra that
are dissimilar from solar produce ratty peaks. We have implemented a
bootstrapping approach to accomodate different spectral types. We
establish a ladder of spectra seperated in effective temperature by
500--1000 K. These spectra are:

::

   $ ls ${SMEMP_WKDIR}/spectra/
   nso.fits # solar spectrum, ~5800 K
   rj55.1872.fits # ~4800 K
   rj26.532.fits # ~3700 K
   rj59.1926.fits # ~6200 K

Shift the spectra by running ``  `` which places
spectra in this directory:

::


   $ python specmatchemp/buildlib/shift_references.py
   Shifting spectrum rj55.1872 onto reference nso
   Shifting spectrum rj26.532 onto reference nso
   Shifting spectrum rj59.1926 onto reference rj55.1872
   
   $ ls ${SMEMP_WKDIR}/shifted_spectra/
   nso_adj.fits
   rj26.532_adj.fits
   rj55.1872_adj.fits
   rj59.1926_adj.fits

Build Table of Library Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The parameters of the stars in the SpecMatch library come from a
variety of literature sources. Run ``python xxx.py`` to read the
parameters into a master csv file. The library sources are
heterogeneous, giving (some have logg while others give stellar
radius). Missing values are given as nans.

The missing values are filled in with Dartmouth isochrones. The once
filled in, the parameters are stored in this file:

::

   ${SMEMP_WKDIR}/xxx.csv


Shift Library Spectra onto Rest Wavelength Scale
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We shift the remainder of the library spectra onto the rest
wavelength scale with the following command.


::

   $ script.py 

Build HDF5 file
~~~~~~~~~~~~~~~

Finally, we assemble all the intermediate files into a monolithic HDF5
file that contains the spectra and their associated parameters using
the following script.

::

   $ script.py
