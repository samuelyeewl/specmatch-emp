#!/usr/bin/env python

"""
Place a library spectrum onto a new wavelength scale with the following properties
- with constant difference in log lambda
- shifted to align with solar (Ganymede) spectrum

@date 05-Apr-2016
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os

def adjust_spectra(path, shift_reference=None):
    """
    Adjusts the given spectrum

    Args:
        path: 
            path to FITS file containing spectrum
        shift_reference:
            path to FITS file containing reference spectrum to shift to
            if None specified, spectrum will not be shifted, only placed onto
            constant log lambda scale.

    Saves the adjusted spectrum to a FITS file with the same name with _adj appended
    e.g. rj76.283_adj.fits
    """

    # open file and read data
    hdu = fits.open(path)
    s = hdu[0].data
    serr = hdu[1].data
    w = hdu[2].data

    # normalize each order to the 95th percentile
    percen_order = np.percentile(s, 95, axis=1)
    s /= percen_order.reshape(-1,1)

    # place spectrum on constant log-lambda wavelength scale
    slog, serrlog, wlog = rescale_w(s, serr, w)

    # solve for velocity shifts between spectra
    if shift_reference is not None:
        # read reference spectrum
        hdu_ref = fits.open(shift_reference)
        s_ref = hdu_ref[0].data
        serr_ref = hdu_ref[1].data
        w_ref = hdu_ref[2].data
        # normalize reference spectrum
        percen_order_ref = np.percentile(s_ref, 95, axis=1)
        s_ref /= percen_order_ref.reshape(-1,1)

        for i in range(len(wlog)):
        # for i in [2]:
            ww = wlog[i]
            ss = slog[i]
            ww_ref = w_ref[i]
            ss_ref = s_ref[i]

            # super sample spectrum to solve for sub-pixel shifts
            logw_min = np.log10(ww[0])
            logw_max = np.log10(ww[-1])
            w_inter = np.logspace(logw_min, logw_max, 10*len(ww), base=10.0)
            # w_inter = np.logspace(logw_min, logw_max, len(ww), base=10.0)
            s_inter = np.interp(w_inter, ww, ss)
            dw = np.median(w_inter[1:] - w_inter[:-1])

            # place reference spectrum on same wavelength scale as library spectrum
            slog_ref = np.interp(w_inter, ww_ref, ss_ref)

            # correlate
            xcorr = np.correlate(slog_ref-1, s_inter-1, mode='same')
            # number of pixels 
            npix = xcorr.shape[0]
            lag_arr = np.arange(-npix/2+1, npix/2+1, 1)
            lag = lag_arr[np.argmax(xcorr)]

            # shift spectrum
            w_inter += lag*dw
            # resample back down to original wavelength spacing
            logw_min = np.log10(w_inter[0])
            logw_max = np.log10(w_inter[-1])
            w_shifted = np.logspace(logw_min, logw_max, len(ww), base=10.0)

            wlog[i] = w_shifted

        plt.plot(wlog[2], slog[2])
        plt.plot(w_ref[2], s_ref[2])
        plt.show()

    # save file
    outfile = os.path.splitext(path)[0] + '_adj.fits'
    hdu[0].data = slog
    hdu[1].data = serrlog
    hdu[2].data = wlog
    hdu.writeto(outfile)

def rescale_w(s, serr, w):
    """
    Place the given spectrum on a constant log-lambda wavelength scale
    """
    wlog = np.empty_like(w)
    slog = np.empty_like(s)
    serrlog = np.empty_like(serr)

    # create a wavelength scale that is uniform in log(lambda)
    for i in range(len(w)):
        ww = w[i]
        ss = s[i]
        sserr = serr[i]

        logw_min = np.log10(ww[0])
        logw_max = np.log10(ww[-1])

        wlog[i] = np.logspace(logw_min, logw_max, len(ww), base=10.0)
        slog[i] = np.interp(wlog[i], ww, ss)
        serrlog[i] = np.interp(wlog[i], ww, sserr)

    return slog, serrlog, wlog

if __name__ == '__main__':
    adjust_spectra('/Users/samuel/Dropbox/SpecMatch-Emp/spectra/iodfitsdb/rj187.477.fits',
        '/Users/samuel/Dropbox/SpecMatch-Emp/spectra/iodfitsdb/rj76.283.fits')