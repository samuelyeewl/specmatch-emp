#!/usr/bin/env python
"""
@filename specmatch_io.py

Provides functions to read and write spectra for use in specmatchemp
"""

import numpy as np
import pandas as pd
from astropy.io import fits
from specmatchemp import pdplus

io_types = ['standard', 'hires']

def read_hires_spectrum(path):
    """
    Opens a spectrum from HIRES.

    Args:
        path: The path to the spectrum
    Returns:
        s: Spectrum
        w: Wavelength scale
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

    return s, w, serr, header

def read_standard_spectrum(path):
    """
    Reads a spectrum stored in the standard format written
    by this module.

    Args:
        path: The path to the spectrum
    Returns:
        s: Spectrum
        w: Wavelength scale
        serr: Spectrum error
        header: File header
    """
    try:
        hdu = fits.open(path)
    except FileNotFoundError:
        print(path+" was not found")
        raise
    
    data = hdu[1].data

    return data['s'], data['w'], data['serr'], hdu[0].header

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



def save_standard_spectrum(path, s, w, serr=None, header=None):
    """
    Saves the given spectrum

    Args:
        path: The path to save the spectrum
        s: Spectrum
        w: Wavelength scale
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
