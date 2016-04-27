#!/usr/bin/env python
"""
Place a target spectrum onto a new wavelength scale onto a new
wavelength scale
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os, sys

def read_spectrum(path):
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


def adjust_spectra(s, serr, w, s_ref, serr_ref, w_ref):
    """
    Adjusts the given spectrum by first placing it on the same wavelength scale as
    the specified reference spectrum, then solves for shifts between the two
    spectra.

    Args:
        s, serr, w: 
            Target spectrum, error and wavelength scale
        s_ref, serr_ref, w_ref:
            Reference spectrum, error, and wavelength scale

    Returns: 
        s_adj, serr_adj, w_adj:
            The adjusted spectrum and wavelength scale.
    """
    # normalize each order of the target spectrum to the 95th percentile
    percen_order = np.percentile(s, 95, axis=1)
    s /= percen_order.reshape(-1,1)

    s_shifted = np.asarray([])
    serr_shifted = np.asarray([])
    ws = np.asarray([])

    # for every order in the spectrum
    for i in range(len(s)):
    # for i in [5]:
        ww = w[i]
        ss = s[i]
        sserr = serr[i]
        # get the target spectrum wavelength range
        w_min = ww[0]
        w_max = ww[-1]
        # remove obvious noise
        mask = np.asarray([True if sp < 1.2 else False for sp in ss])
        ss = ss[mask]
        sserr = sserr[mask]
        ww = ww[mask]

        # get the reference spectrum in the same range
        in_range = np.asarray([True if wr > w_min and wr < w_max else False
            for wr in w_ref])
        start_idx = np.argmax(in_range)
        w_ref_c = w_ref[in_range]
        s_ref_c = s_ref[in_range]

        # place the target spectrum on the same wavelength scale
        ss, sserr = rescale_w(ss, sserr, ww, w_ref_c)

        # solve for shifts in different sections
        num_sections = 8
        l_sect = int(len(w_ref_c)/num_sections)      # length of each section
        lags = np.empty(num_sections)
        center_pix = np.empty(num_sections)

        for j in range(num_sections):
            # get the shifts in pixel number
            lag, lag_arr, xcorr = solve_for_shifts(ss[j*l_sect:(j+1)*l_sect],
                s_ref_c[j*l_sect:(j+1)*l_sect])
            lags[j] = lag
            center_pix[j] = (j+1/2)*l_sect
            # plt.plot(lag_arr, xcorr)

        # we expect that the shifts across every order are close together
        # so we should remove outliers
        med = np.median(lags)
        tol = 1     # permitted deviation from median
        not_outlier = np.asarray([True if l > med-tol and l < med+tol else False 
            for l in lags])
        lags = lags[not_outlier]
        center_pix = center_pix[not_outlier]

        # fit a straight line to the shifts
        fit = np.polyfit(center_pix, lags, 1)
        pix_arr = np.arange(0, len(w_ref_c))
        pix_shifted = pix_arr - fit[1] - pix_arr*fit[0]

        # plt.plot(center_pix, lags)
        # plt.plot(pixs, fit[0]*pixs+fit[1])
        # plt.show()

        # plt.plot(pixs, pix_shifted)
        # plt.show()

        # now shift the spectrum
        # ww_shifted = w_ref_c - fit[1] - w_ref_c*fit[0]
        # plt.plot(ww_shifted, w_ref_c)
        # plt.show()

        pix_min = int(pix_shifted[0]) + 1
        pix_max = int(pix_shifted[-1])
        # new pixel array
        new_pix = np.arange(pix_min, pix_max)
        # new wavelength array
        w_ref_c = w_ref[start_idx+pix_min:start_idx+pix_max]

        # interpolate the spectrum back onto the reference spectrum
        ss_shifted = np.interp(new_pix, pix_shifted, ss)
        sserr_shifted = np.interp(new_pix, pix_shifted, sserr)

        # append to array, discarding leading and trailing 50 pixels
        s_shifted = np.append(s_shifted, ss_shifted[50:-50])
        serr_shifted = np.append(serr_shifted, sserr_shifted[50:-50])
        ws = np.append(ws, w_ref_c[50:-50])

    # average any values that appear twice
    w_flattened = np.unique(ws)
    s_flattened = np.empty_like(w_flattened)
    serr_flattened = np.empty_like(w_flattened)

    for i, wl in enumerate(w_flattened):
        s_flattened[i] = np.mean(s_shifted[ws == wl])
        serr_flattened[i] = np.mean(serr_shifted[ws == wl])

    return s_flattened, serr_flattened, w_flattened
    
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
    snew = np.interp(w_ref, w, s)
    serrnew = np.interp(w_ref, w, s)

    # snew = np.empty_like(s)
    # serrnew = np.empty_like(serr)

    # for i in range(len(w)):
    #     snew[i] = np.interp(w_ref[i], w[i], s[i])
    #     serrnew[i] = np.interp(w_ref[i], w[i], serr[i])

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
