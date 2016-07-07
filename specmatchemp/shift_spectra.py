#!/usr/bin/env python
"""
@filename shift_spectra.py

Shift a target spectrum onto a reference spectrum.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import argparse
from specmatchemp.io import specmatchio

def adjust_spectra(s, serr, w, s_ref, serr_ref, w_ref, diagnostic=False, outfile='./diag.csv', diagnostic_hdr=None):
    """
    Adjusts the given spectrum by first placing it on the same wavelength scale as
    the specified reference spectrum, then solves for shifts between the two
    spectra.

    Args:
        s, serr, w: 
            Target spectrum, error and wavelength scale
        s_ref, serr_ref, w_ref:
            Reference spectrum, error, and wavelength scale
        diagnostic:
            Set to true if diagnostic plots of lags are required

    Returns: 
        s_adj, serr_adj, w_adj:
            The adjusted spectrum and wavelength scale.
    """
    # normalize each order of the target spectrum by fitting a spline
    percen_order = np.percentile(s, 95, axis=1)
    s /= percen_order.reshape(-1,1)
    serr /= percen_order.reshape(-1,1)

    s_shifted = np.asarray([[]])
    serr_shifted = np.asarray([[]])
    ws = np.asarray([[]])

    tol = 5     # permitted deviation from median
    num_sections = 7    # number of sections of each order

    if diagnostic:
        f = open(outfile, 'w')
        if diagnostic_hdr is not None:
            f.write(diagnostic_hdr)
            f.write('# tol = {0:d}, num_sections = {1:d}\n'.format(tol, num_sections))
        f.close()
        f = open(outfile, 'ab')
        list_lags=np.empty((len(s),3,num_sections))

    # for every order in the spectrum
    for i in range(len(s)):
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
        l_sect = int(len(w_ref_c)/num_sections)
        lags = np.empty(num_sections)
        center_pix = np.empty(num_sections)

        for j in range(num_sections):
            # get the shifts in pixel number
            lag, lag_arr, xcorr = solve_for_shifts(ss[j*l_sect:(j+1)*l_sect],
                s_ref_c[j*l_sect:(j+1)*l_sect])
            lags[j] = lag
            center_pix[j] = (j+1/2)*l_sect

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

        # if(diagnostic):
        #     plt.plot(center_pix, lags)
        #     plt.plot(pix_arr, fit[0]*pix_arr+fit[1])

        if diagnostic:
            fitted = fit[0]*center_pix+fit[1]
            list_lags[i] = [center_pix, lags, fitted]

        # plt.plot(pixs, pix_shifted)
        # plt.show()

        # now shift the spectrum
        # ww_shifted = w_ref_c - fit[1] - w_ref_c*fit[0]
        # plt.plot(ww_shifted, w_ref_c)
        # plt.show()

        pix_min = max(int(pix_shifted[0]) + 1, 0)
        pix_max = min(int(pix_shifted[-1]), len(w_ref_c))

        # new pixel array
        new_pix = np.arange(pix_min, pix_max)
        # new wavelength array
        w_ref_c = w_ref[start_idx+pix_min:start_idx+pix_max]
        
        # interpolate the spectrum back onto the reference spectrum
        ss_shifted = np.interp(new_pix, pix_shifted, ss)
        sserr_shifted = np.interp(new_pix, pix_shifted, sserr)

        # append to array, discarding leading and trailing 20 pixels
        s_shifted = np.append(s_shifted, ss_shifted[20:-20])
        serr_shifted = np.append(serr_shifted, sserr_shifted[20:-20])
        ws = np.append(ws, w_ref_c[20:-20])

    if diagnostic:
        # reshape into 2d array for storage
        list_lags = list_lags.reshape((len(s),3*num_sections))
        np.savetxt(f, list_lags)
        f.close()

    return s_shifted, serr_shifted, ws

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
    xcorr = np.correlate(s-1, s_ref-1, mode='same')
    xcorr = np.nan_to_num(xcorr)
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
    serrnew = np.interp(w_ref, w, serr)

    return snew, serrnew

def shift(target_path, target_type, reference_path, output_path, diagnostic=False, diagnosticfile='img.img'):
    if target_type =='hires':
        try:
            s, w, serr, header = specmatch_io.read_hires_spectrum(target_path)
        except:
            raise
    elif target_type =='standard':
        try:
            s, w, serr, header = specmatch_io.read_standard_spectrum(target_path)
        except:
            raise

        s = np.asarray([s])
        w = np.asarray([w])
        serr = np.asarray([serr])

    try:
        s_ref, w_ref, serr_ref, header_ref = specmatch_io.read_standard_spectrum(reference_path)
    except:
        raise

    s_adj, serr_adj, w_adj = adjust_spectra(s, serr, w, s_ref, serr_ref, w_ref, diagnostic, diagnosticfile)

    try:
        specmatch_io.save_standard_spectrum(output_path, s_adj, w_adj, serr_adj, header)
    except:
        print('Could not save to '+output_path)
        raise

    return s_adj, w_adj, serr_adj


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Shift a target spectrum onto a reference spectrum')
    parser.add_argument('target_path', type=str)
    parser.add_argument('target_type', choices=io_types)
    parser.add_argument('reference_path', type=str)
    parser.add_argument('output_path', type=str)
    args = parser.parse_args()

    shift(args.target_path, args.target_type, args.reference_path, args.output_path)