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
~XXX K. These spectra are:

::

    ${SMEMP_WKDIR}/spectra/
        rjXXX.fits
	rjXXX.fits


Shift the spectra by running ``python create_wavref.py`` which places
spectra in this directory:

::

    ${SMEMP_WKDIR}/spectra_restwav/
        rjXXX_restwav.fits
	rjXXX_restwav.fits


Build Table of Library Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The parameters of the stars in the SpecMatch library come from a
variety of literature sources. Run ``python xxx.py`` to read the
parameters into a master csv file. The library sources are
heterogeneous, giving (some have logg while others give stellar
radius). Missing values are given as nans.

The missing values are filled in with Dartmouth isochrones. The once filled in, the parameters are stored in this file:

::

   ${SMEMP_WKDIR}/






