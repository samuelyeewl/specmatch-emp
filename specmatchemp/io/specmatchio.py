"""
@filename specmatch.io.specmatch_io.py
Provides functions to read and write spectra for use in specmatchemp
"""

from __future__ import print_function

import numpy as np
import pandas as pd
from astropy.io import fits
from specmatchemp.io import pdplus
from specmatchemp.spectrum import Spectrum

io_types = ['standard', 'hires']

def read_hires_spectrum(path):
    """
    Opens a spectrum from HIRES.
    Args:
        path: The path to the spectrum
    Returns:
        w: Wavelength scale
        s: Spectrum
        serr: Spectrum error
        header: File header
    """
    try:
        hdu = fits.open(path)
    except FileNotFoundError:
        print(path+" was not found")
        raise

    s = hdu[0].data
    serr = hdu[1].data
    w = hdu[2].data
    header = hdu[0].header

    return w, s, serr, header

def read_standard_spectrum(path, wavlim=None):
    """
    Reads a spectrum stored in the standard format written
    by this module.
    Args:
        path: The path to the spectrum
        wavlim (2-element iterable): (optional) Fix the wavelength range
            to be read.
    Returns:
        w: Wavelength scale
        s: Spectrum
        serr: Spectrum error
        header: File header
    """
    try:
        hdu = fits.open(path)
    except FileNotFoundError:
        print(path+" was not found")
        raise
    
    data = hdu[1].data
    w = data['w']
    s = data['s']
    serr = data['serr']

    if wavlim is not None:
        wavidx, = np.where((w > wavlim[0]) & (w < wavlim[1]))
        idxmin = wavidx[0]
        idxmax = wavidx[-1]+1
        w = w[idxmin:idxmax]
        s = s[idxmin:idxmax]
        serr = s[idxmin:idxmax]

    return w, s, serr, hdu[0].header

def read_as_dataframe(path):
    """
    Reads a spectrum stored in the standard format written by
    this module into a dataframe.
    Args:
        path: The path to the spectrum
    Returns:
        df: Dataframe containing the data
        header: File header
    """

    try:
        hdu = fits.open(path)
    except FileNotFoundError:
        print(path+" was not found")
        raise

    data = hdu[1].data
    df = pd.DataFrame(dict(s=data['s'],w=data['w'], serr=data['serr']))
    df = pdplus.LittleEndian(df.to_records(index=False))
    df = pd.DataFrame(df)

    return df, hdu[0].header

def truncate_spectrum(wavlim, w, s, serr=None):
    """Truncates a spectrum to be within the given wavelength limits.

    Args:
        wavlim (2-element iterable): The wavelength limits
        w (np.ndarray): Wavelength array
        s (np.ndarray) : Spectrum
        serr (np.ndarray) : (optional) Uncertainty in spectrum
    Returns:
        w_truncated, s_truncated, serr_truncated
    """
    wavidx, = np.where((w >= wavlim[0]) & (w <= wavlim[1]))
    idxmin = wavidx[0]
    idxmax = wavidx[-1]+1
    w = w[idxmin:idxmax]
    s = s[idxmin:idxmax]
    if serr is not None:
        serr = serr[idxmin:idxmax]
        return w, s, serr
    else:
        return w, s

def save_standard_spectrum(path, w, s, serr=None, header=None):
    """
    Saves the given spectrum
    Args:
        path: The path to save the spectrum
        w: Wavelength scale
        s: Spectrum
        serr: Spectrum error
        header: File header (optional)
    """
    # Create a header
    if header is None:
        header = fits.Header()
    prihdu = fits.PrimaryHDU(header=header)

    # If no error is provided
    if serr is None:
        serr = np.empty_like(s)
        serr.fill(None)

    tbhdu = fits.BinTableHDU.from_columns(
        [fits.Column(name='s', format='D', array=s),
        fits.Column(name='w', format='D', array=w),
        fits.Column(name='serr', format='D', array=serr)])

    hdulist = fits.HDUList([prihdu, tbhdu])
    hdulist.writeto(path, clobber=True)