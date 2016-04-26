#!/usr/bin/env python
"""Place a target spectrum onto a new wavelength scale onto a new
wavelength scale
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os, sys

def open_spectrum(path):
    """
    Opens the spectrum at the specified path.

    Returns:
        s, serr, w
    """
    # open file and read data
    hdu = fits.open(path)
    s = hdu[0].data
    serr = hdu[1].data
    w = hdu[2].data

    return s, serr, w, hdu


def adjust_spectra(s, serr, w, s_ref=None, serr_ref=None, w_ref=None):
    """
    Adjusts the given spectrum by placing it on a constant log-lambda scale.
    If a reference spectrum is specified, solves for shifts between target
    and reference spectra and shifts the target spectrum accordingly.

    Args:
        s, serr, w: 
            Target spectrum, error and wavelength scale
        s_ref, serr_ref, w_ref:
            Reference spectrum, error, and wavelength scale (optional)

    Returns: 
        s_adj, serr_adj, w_adj:
            The adjusted spectrum and wavelength scale.
    """
    # normalize each order of the target spectrum to the 95th percentile
    percen_order = np.percentile(s, 95, axis=1)
    s /= percen_order.reshape(-1,1)

    # place spectrum on constant log-lambda wavelength scale
    slog, serrlog, wlog = rescale_log_w(s, serr, w)

    # solve for velocity shifts between spectra
    if s_ref is not None:
        w_shifted = np.empty_like(wlog)

        # normalize reference spectrum
        percen_order_ref = np.percentile(s_ref, 95, axis=1)
        s_ref /= percen_order_ref.reshape(-1,1)

        # place reference spectrum on same wavelength scale as target spectrum
        s_ref, serr_ref = rescale_w(s_ref, serr_ref, w_ref, wlog)

        # for i in range(len(wlog)):
        for i in [2]:
            ww = wlog[i]
            ss = slog[i]
            ww_ref = w_ref[i]
            ss_ref = s_ref[i]

            # solve for shifts in different sections
            num_sections = 2
            lags = np.empty(num_sections)
            l_sect = len(ss)/num_sections
            ww_shifted = np.empty_like(ww)

            for j in range(num_sections):
                # split array
                ss_sect = ss[j*l_sect:(j+1)*l_sect]
                ss_ref_sect = ss_ref[j*l_sect:(j+1)*l_sect]
                
                # solve for shifts
                lag, lag_arr, xcorr = solve_for_shifts(ss_sect, ss_ref_sect)

                # shift spectrum
                dw = np.median(ww[j*l_sect+1:(j+1)*l_sect] - ww[j*l_sect:(j+1)*l_sect-1])
                ww_shifted[j*l_sect:(j+1)*l_sect] = ww[j*l_sect:(j+1)*l_sect] - dw*lag

                lags[j] = lag

                plt.plot(lag_arr, xcorr)
            
            plt.show()

            # plot variation of lags
            plt.plot(np.linspace(ww[0], ww[-1], num_sections), lags)
            plt.show()

            w_shifted[i] = ww_shifted

    else:
        w_shifted = wlog

    return slog, serrlog, w_shifted

def solve_for_shifts(s, s_ref):
    """
    Solve for the pixel shifts required to align two spectra that are on the same
    wavelength scale.

    Correlates the two spectra, then fits a quadratic to the peak in order to
    solve for sub-pixel shifts.

    Args:
        s: The target spectrum
        s_ref: The reference spectrum
        w: The common wavelength scale.
    
    Returns:
        The pixel shift, the lag and correlation data
    """
    # correlate the two spectra
    xcorr = np.correlate(s-1, s_ref-1, mode='same')
    max_corr = np.argmax(xcorr)

    # number of pixels
    npix = xcorr.shape[0]
    lag_arr = np.arange(-npix/2+1, npix/2+1, 1)

    # select points around the peak and fit a quadratic
    lag_peaks = lag_arr[max_corr-5:max_corr+5]
    xcorr_peaks = xcorr[max_corr-5:max_corr+5]
    p = np.polyfit(lag_peaks, xcorr_peaks, 2)
    # peak is simply -p[1]/2p[0]
    lag = -p[1]/(2*p[0])

    return lag, lag_arr, xcorr


def rescale_w(s, serr, w, w_ref):
    """
    Place the given spectrum on the wavelength scale specified by w_ref

    Args:
        s, serr, w: The spectrum and original wavelength scale.
        w_ref: The desired wavelength scale

    Returns:
        The spectrum and associated error on the desired scale.
    """

    snew = np.empty_like(s)
    serrnew = np.empty_like(serr)

    for i in range(len(w)):
        snew[i] = np.interp(w_ref[i], w[i], s[i])
        serrnew[i] = np.interp(w_ref[i], w[i], serr[i])

    return snew, serrnew


def rescale_log_w(s, serr, w):
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
