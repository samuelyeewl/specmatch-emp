#!/usr/bin/env python
"""
@filename shift_spectra.py

Shift a target spectrum onto a reference spectrum.
"""
import numpy as np
import matplotlib.pyplot as plt
import argparse
from specmatchemp.io import specmatchio

def adjust_spectra(s, serr, w, s_ref, serr_ref, w_ref, diagnostic=False, diagnosticfile='img.img'):
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
    # normalize each order of the target spectrum to the 95th percentile
    percen_order = np.percentile(s, 95, axis=1)
    s /= percen_order.reshape(-1,1)

    s_shifted = np.asarray([])
    serr_shifted = np.asarray([])
    ws = np.asarray([])

    plt.cla()

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

        # solve for shifts in different sections (of length 1200 pix)
        l_sect = 1000
        num_sections = int(len(w_ref_c)/l_sect)
        l_sect = int(len(w_ref_c)/num_sections)
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
        tol = 2     # permitted deviation from median
        not_outlier = np.asarray([True if l > med-tol and l < med+tol else False 
            for l in lags])
        lags = lags[not_outlier]
        center_pix = center_pix[not_outlier]

        # fit a straight line to the shifts
        fit = np.polyfit(center_pix, lags, 1)
        pix_arr = np.arange(0, len(w_ref_c))
        pix_shifted = pix_arr - fit[1] - pix_arr*fit[0]

        # if (diagnostic and i==2):
        if(diagnostic):
            plt.plot(center_pix, lags)
            plt.plot(pix_arr, fit[0]*pix_arr+fit[1])
        # plt.show()

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

        # append to array, discarding leading and trailing 50 pixels
        s_shifted = np.append(s_shifted, ss_shifted[50:-50])
        serr_shifted = np.append(serr_shifted, sserr_shifted[50:-50])
        ws = np.append(ws, w_ref_c[50:-50])

    if(diagnostic):
        plt.savefig(diagnosticfile)
        plt.clf()

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