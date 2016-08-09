#!/usr/bin/env python
"""
@filename shift.py

Shift a target spectrum onto a reference spectrum.
"""
import argparse

import h5py
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.optimize import least_squares

from specmatchemp.spectrum import Spectrum
from specmatchemp.io import specmatchio

def shift(targ, ref, mask=None, outfile=None):
    """Shifts the given spectrum by placing it on the same wavelength
    scale as the specified reference spectrum, then solves for shifts
    between the two spectra through cross-correlation.

    Args:
        targ (Spectrum): Target spectrum
        ref (Spectrum): Reference spectrum
        mask (optional [pd.DataFrame]): Mask for telluric lines
        outfile (optional [file]): h5 file to store diagnostic data in

    Returns:
        shifted (Spectrum): Adjusted and flattened spectrum
    """
    s = targ.s
    serr = targ.serr
    w = targ.w

    # normalize each order of the target spectrum by dividing by the 95th percentile
    percen_order = np.percentile(s, 95, axis=1)
    s /= percen_order.reshape(-1,1)

    # create empty 2d arrays to store each order
    s_shifted = np.asarray([[]])
    serr_shifted = np.asarray([[]])
    mask_shifted = np.asarray([[]])
    ws = np.asarray([[]])

    # create lists to store diagnostic data
    lag_data = []
    center_pix_data = []
    fit_data = []

    num_sections = 7

    # shift each order
    for i in range(s.shape[0]):
        ww = w[i]
        ss = s[i]
        sserr = serr[i]

        # create mask
        mm = np.empty_like(ww)
        mm.fill(True)
        if mask is not None:
            mask_cut = mask.query('order == {0:d}'.format(i))
            for n, row in mask_cut.iterrows():
                start = row['minpix']
                end = row['maxpix']
                mm[start:end] = False

        if i == 15:
            import matplotlib.pyplot as plt
            plt.plot(ss-2)
            plt.plot(mm-2)


        # clip ends off each order
        cliplen = 15
        ww = ww[cliplen:-cliplen]
        ss = ss[cliplen:-cliplen]
        sserr = sserr[cliplen:-cliplen]
        mm = mm[cliplen:-cliplen]

        # clip obvious noise
        clip = np.asarray([True if sp < 1.2 else False for sp in ss])
        ss = ss[clip]
        sserr = sserr[clip]
        mm = mm[clip]
        ww = ww[clip]

        # get the reference spectrum in the same range as the target range
        w_min = ww[0]
        w_max = ww[-1]

        in_range = np.asarray([True if wr > w_min and wr < w_max else False
            for wr in ref.w])
        start_idx = np.argmax(in_range)
        w_ref_c = ref.w[in_range]
        s_ref_c = ref.s[in_range]

        # place the target spectrum on the same wavelength scale
        ss, sserr, mm = rescale_w(ss, sserr, ww, mm, w_ref_c)

        if i == 15:
            plt.plot(ss-1)
            plt.plot(mm-1)

        # solve for shifts in different sections
        l_sect = int(len(w_ref_c)/num_sections)
        lags = np.empty(num_sections)
        center_pix = np.empty(num_sections)
        lag_arrs = []
        xcorrs = []

        if outfile is not None:
            grp = outfile.create_group("/xcorr/order_{0:d}".format(i))

        for j in range(num_sections):
            ss_sect = ss[j*l_sect:(j+1)*l_sect]
            s_ref_sect = s_ref_c[j*l_sect:(j+1)*l_sect]
            # get the shifts in pixel number
            lag, lag_arr, xcorr = solve_for_shifts(ss_sect, s_ref_sect)
            lags[j] = lag
            center_pix[j] = (j+1/2)*l_sect

            if outfile is not None:
                subgrp = grp.create_group("sect_{0:d}".format(j))
                subgrp["xcorr"] = xcorr
                subgrp["lag_arr"] = lag_arr

        # remove clear outliers
        med = np.median(lags)
        tol = 20
        not_outlier = np.asarray([True if l > med-tol and l < med+tol else False for l in lags])
        lags_trunc = lags[not_outlier]
        center_pix_trunc = center_pix[not_outlier]

        # use robust least squares to fit a line to the shifts (Cauchy loss function)
        p_guess = np.array([0,0])
        fit_res = least_squares(_linear_fit_residuals, p_guess, \
            args=(center_pix_trunc, lags_trunc), loss='cauchy')
        fit = fit_res.x
        pix_arr = np.arange(0, len(w_ref_c))
        pix_shifted = pix_arr - fit[1] - pix_arr*fit[0]

        # don't read past the wavelength array
        pix_min = max(int(pix_shifted[0]), 0)
        pix_max = min(int(pix_shifted[-1]), len(ref.w)-start_idx)

        # new pixel array
        new_pix = np.arange(pix_min, pix_max)
        # new wavelength array
        w_ref_c = ref.w[start_idx+pix_min:start_idx+pix_max]

        # interpolate the spectrum back onto the reference spectrum
        ss_shifted = np.interp(new_pix, pix_shifted, ss)
        sserr_shifted = np.interp(new_pix, pix_shifted, sserr)
        mm_shifted = np.interp(new_pix, pix_shifted, mm)


        # append to array
        s_shifted = np.append(s_shifted, ss_shifted)
        serr_shifted = np.append(serr_shifted, sserr_shifted)
        mask_shifted = np.append(mask_shifted, mm_shifted)
        ws = np.append(ws, w_ref_c)

        if i == 15:
            plt.plot(ss)
            # shifted mask
            plt.plot(ss_shifted+1)
            plt.plot(mm_shifted+1)

        # save diagnostic data
        lag_data.append(lags)
        center_pix_data.append(center_pix)
        fitted = fit[0]*center_pix+fit[1]
        fit_data.append(fitted)


    # save diagnostic data
    if outfile is not None:
        outfile.create_dataset('lag', data=np.asarray(lag_data))
        outfile.create_dataset('center_pix', data=np.asarray(center_pix_data))
        outfile.create_dataset('fit', data=np.asarray(fit_data))

    # flatten spectrum
    w_min = ws[0]
    w_max = ws[-1]
    in_range = np.asarray([True if wr > w_min and wr < w_max else False for wr in ref.w])
    w_ref_trunc = ref.w[in_range]

    w_flat, s_flat, serr_flat, mask_flat = flatten(ws, s_shifted, serr_shifted, mask_shifted, w_ref=w_ref_trunc)

    return Spectrum(w_flat, s_flat, serr_flat, name=targ.name, mask=mask_flat, header=targ.header, attrs=targ.attrs)


def _isclose(a, b, abs_tol=1e-6):
    """Small helper function to determine if two floats are close.
    Only accepts absolute tolerances.
    """
    return abs(a-b) <= abs_tol

def _fill_nans(s, fill):
    """Replaces nans with the provided fill value
    """
    s[np.isnan(s)] = fill

    return s

def _linear_fit_residuals(p, x, y):
    """Calculates residuals for a linear fit
    """
    return p[0]*x + p[1] - y


def flatten(w, s, serr, mask, w_ref=None, wavlim=None):
    """Flattens a given 2-D spectrum into a 1-D array.
    Merges overlapped points by taking the mean.
    If w_ref is given, fills values that don't occur in the 2D spectrum
    with np.nan

    Args:
        w (np.ndarray): Wavelength array
        s (np.ndarray): Spectrum
        serr (np.ndarray): Uncertainty in spectrum
        mask (np.ndarray): Boolean mask
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
    mask_flattened = np.empty_like(w_flattened)
    idx_max = len(w)-1
    c_idx = 0
    n_idx = 0

    for i, wl in enumerate(w_ref):
        while w[c_idx] < wl and c_idx < idx_max and not _isclose(w[c_idx], wl):
            c_idx += 1

        if c_idx >= idx_max:
            s_flattened[i] = np.nan
            serr_flattened[i] = np.nan
            mask_flattened[i] = False
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
                mask_flattened[i] = bool(mask[c_idx]) & bool(mask[n_idx])
            else:
                s_flattened[i] = s[c_idx]
                serr_flattened[i] = serr[c_idx]
                mask_flattened[i] = mask[c_idx]
            c_idx += 1
        else:
            s_flattened[i] = np.nan
            serr_flattened[i] = np.nan
            mask_flattened[i] = False

    return w_flattened, s_flattened, serr_flattened, mask_flattened

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
    smean = np.nanmean(s)
    srefmean = np.nanmean(s_ref)
    # fill nans with mean value so they contribute nothing to correlation
    s = _fill_nans(s, smean)
    s_ref = _fill_nans(s_ref, srefmean)

    xcorr = np.correlate(s-smean, s_ref-srefmean, mode='full')
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

def rescale_w(s, serr, w, m, w_ref):
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
    mnew = np.interp(w_ref, w, m).astype(bool)

    return snew, serrnew, mnew


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Shift a target spectrum onto a reference spectrum')
    parser.add_argument('target_path', type=str)
    parser.add_argument('target_type', choices=io_types)
    parser.add_argument('reference_path', type=str)
    parser.add_argument('output_path', type=str)
    args = parser.parse_args()

    shift(args.target_path, args.target_type, args.reference_path, args.output_path)