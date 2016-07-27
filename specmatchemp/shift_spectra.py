#!/usr/bin/env python
"""
@filename shift_spectra.py

Shift a target spectrum onto a reference spectrum.
"""
import argparse

import h5py
import numpy as np
from scipy.interpolate import UnivariateSpline

from specmatchemp.io import specmatchio

def shift(s, serr, w, s_ref, serr_ref, w_ref, outfile=None):
    """Shifts the given spectrum by placing it on the same wavelength
    scale as the specified reference spectrum, then solves for shifts
    between the two spectra through cross-correlation.

    Args:
        s, serr, w : Target spectrum, error and wavelength scale
        s_ref, serr_ref, w_ref : Reference spectrum, error and wavelength scale
        outfile : (optional) h5 file to store diagnostic data in

    Returns:
        s_adj, serr_adj, w_adj: Adjusted and flattened spectrum
    """
    # normalize each order of the target spectrum by dividing by the 95th percentile
    percen_order = np.percentile(s, 95, axis=1)
    s /= percen_order.reshape(-1,1)

    # create empty 2d arrays to store each order
    s_shifted = np.asarray([[]])
    serr_shifted = np.asarray([[]])
    ws = np.asarray([[]])

    # create lists to store diagnostic data
    lag_data = []
    center_pix_data = []
    fit_data = []

    tol = 5
    num_sections = 7

    # shift each order
    for i in range(s.shape[0]):
        ww = w[i]
        ss = s[i]
        sserr = serr[i]

        # clip ends off each order
        cliplen = 15
        ww = ww[cliplen:-cliplen]
        ss = ss[cliplen:-cliplen]
        sserr = sserr[cliplen:-cliplen]

        # clip obvious noise
        mask = np.asarray([True if sp < 1.2 else False for sp in ss])
        ss = ss[mask]
        sserr = sserr[mask]
        ww = ww[mask]

        # get the reference spectrum in the same range as the target range
        w_min = ww[0]
        w_max = ww[-1]

        in_range = np.asarray([True if wr > w_min and wr < w_max else False
            for wr in w_ref])
        start_idx = np.argmax(in_range)
        w_ref_c = w_ref[in_range]
        s_ref_c = s_ref[in_range]

        # place the target spectrum on the same wavelength scale
        ss, sserr = rescale_w(ss, sserr, ww, w_ref_c)

        # solve for shifts in different sections
        l_sect = int(len(w_ref_c)/num_sections)
        lags = np.empty(num_sections)
        center_pix = np.empty(num_sections)
        lag_arrs = []
        xcorrs = []

        for j in range(num_sections):
            # get the shifts in pixel number
            lag, lag_arr, xcorr = solve_for_shifts(ss[j*l_sect:(j+1)*l_sect],
                s_ref_c[j*l_sect:(j+1)*l_sect])
            lags[j] = lag
            center_pix[j] = (j+1/2)*l_sect

        lag_data.append(lags)
        center_pix_data.append(center_pix)

        # we expect that the shifts across every order are close together
        # so we should remove outliers
        med = np.median(lags)
        not_outlier = np.asarray([True if l > med-tol and l < med+tol else False 
            for l in lags])
        lags_trunc = lags[not_outlier]
        center_pix_trunc = center_pix[not_outlier]

        # fit a straight line to the shifts
        fit = np.polyfit(center_pix_trunc, lags_trunc, 1)
        pix_arr = np.arange(0, len(w_ref_c))
        pix_shifted = pix_arr - fit[1] - pix_arr*fit[0]

        fitted = fit[0]*center_pix+fit[1]
        fit_data.append(fitted)

        # don't read past the wavelength array
        pix_min = max(int(pix_shifted[0]), 0)
        pix_max = min(int(pix_shifted[-1]), len(w_ref)-start_idx)

        # new pixel array
        new_pix = np.arange(pix_min, pix_max)
        # new wavelength array
        w_ref_c = w_ref[start_idx+pix_min:start_idx+pix_max]

        # interpolate the spectrum back onto the reference spectrum
        ss_shifted = np.interp(new_pix, pix_shifted, ss)
        sserr_shifted = np.interp(new_pix, pix_shifted, sserr)

        # append to array
        s_shifted = np.append(s_shifted, ss_shifted)
        serr_shifted = np.append(serr_shifted, sserr_shifted)
        ws = np.append(ws, w_ref_c)

    # save diagnostic data
    if outfile is not None:
        outfile.create_dataset('lag', data=np.asarray(lag_data))
        outfile.create_dataset('center_pix', data=np.asarray(center_pix_data))
        outfile.create_dataset('fit', data=np.asarray(fit_data))

    # flatten spectrum
    w_min = ws[0]
    w_max = ws[-1]
    in_range = np.asarray([True if wr > w_min and wr < w_max else False for wr in w_ref])
    w_ref_trunc = w_ref[in_range]

    w_flat, s_flat, serr_flat = flatten(ws, s_shifted, serr_shifted, w_ref=w_ref_trunc)

    return s_flat, serr_flat, w_flat


def _isclose(a, b, abs_tol=1e-6):
    """Small helper function to determine if two floats are close.
    Only accepts absolute tolerances.
    """
    return abs(a-b) <= abs_tol

def flatten(w, s, serr, w_ref=None, wavlim=None):
    """Flattens a given 2-D spectrum into a 1-D array.
    Merges overlapped points by taking the mean.
    If w_ref is given, fills values that don't occur in the 2D spectrum
    with np.nan

    Args:
        w (np.ndarray): Wavelength array
        s (np.ndarray): Spectrum
        serr (np.ndarray): (optional) Uncertainty in spectrum
        w_ref (np.nadarray): (optional) Reference, 1-D wavelength array
        wavlim (2-element iterable): (optional) Wavelength limits
    
    Returns:
        w, s, serr: Wavelength, spectrum and uncertainty in spectrum
    """
    assert np.shape(w) == np.shape(s), "w, s not the same shape!"
    assert np.shape(w) == np.shape(serr), "w, serr not the same shape!"

    if w_ref is None:
        w_flattened = np.unique(w)
    else:
        w_flattened = w_ref

    # create new arrays to contain spectrum
    s_flattened = np.empty_like(w_flattened)
    serr_flattened = np.empty_like(w_flattened)
    idx_max = len(w)-1
    c_idx = 0
    n_idx = 0

    for i, wl in enumerate(w_ref):
        while w[c_idx] < wl and c_idx < idx_max and not _isclose(w[c_idx], wl):
            c_idx += 1

        if c_idx >= idx_max:
            s_flattened[i] = np.nan
            serr_flattened[i] = np.nan
            continue

        overlap = False
        if _isclose(w[c_idx], wl):
            # scan for overlapping region
            while n_idx < idx_max:
                n_idx += 1
                if c_idx == n_idx:
                    continue
                elif _isclose(w[n_idx], wl):
                    overlap=True
                    break
                elif w[n_idx] > w[c_idx] and w[n_idx] < w[n_idx-1]:
                    # lock at start of new order
                    n_idx -= 1
                    break
            if overlap:
                s_flattened[i] = (s[c_idx]+s[n_idx])/2
                serr_flattened[i] = (serr[c_idx]+serr[n_idx])/2
            else:
                s_flattened[i] = s[c_idx]
                serr_flattened[i] = serr[c_idx]
            c_idx += 1
        else:
            s_flattened[i] = np.nan
            serr_flattened[i] = np.nan

    return w_flattened, s_flattened, serr_flattened

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
    xcorr = np.correlate(s-1, s_ref-1, mode='full')
    xcorr = np.nan_to_num(xcorr)
    max_corr = np.argmax(xcorr)

    # number of pixels
    npix = xcorr.shape[0]
    lag_arr = np.arange(-(npix-1)/2, (npix+1)/2, 1)

    # select points around the peak and fit a quadratic
    lag_peaks = lag_arr[max_corr-5:max_corr+6]
    xcorr_peaks = xcorr[max_corr-5:max_corr+6]

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
    serrnew = np.interp(w_ref, w, serr)

    return snew, serrnew


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Shift a target spectrum onto a reference spectrum')
    parser.add_argument('target_path', type=str)
    parser.add_argument('target_type', choices=io_types)
    parser.add_argument('reference_path', type=str)
    parser.add_argument('output_path', type=str)
    args = parser.parse_args()

    shift(args.target_path, args.target_type, args.reference_path, args.output_path)