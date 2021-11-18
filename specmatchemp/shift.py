"""
@filename shift.py

Shift a target spectrum onto a reference spectrum.
"""
from __future__ import print_function

import numpy as np
from astropy.io import fits
from scipy.optimize import least_squares
from scipy.special import expit
from scipy.stats import sigmaclip

from specmatchemp import spectrum
from specmatchemp.utils import utils


def bootstrap_shift(targ, ref_list, store=None):
    """Shift a target spectrum using a bootstrapping approach.

    First performs a cross-correlation of the target spectrum in a fixed
    wavelength region against each of the provided references. The reference
    spectrum which gives the highest median cross-correlation peak is used
    to shift the rest of the target spectrum.

    Args:
        targ (spectrum.Spectrum): Target spectrum
        ref_list (list of spectrum.Spectrum): List of reference spectra
        store (dict-like object, optional): h5 file or dict to record
            diagnostic data.

    Returns:
        shifted (Spectrum): Shifted and flattened spectrum.
    """
    print("Shifting spectrum {0}".format(targ.name))
    # Try to use the Mg triplet to determine which reference spectrum is
    # best.
    if isinstance(targ, spectrum.HiresSpectrum) and targ.w.ndim == 2:
        # In the HIRES echelle object, the Mgb triplet is in order 2
        ref_order = 2
        targ_cut = spectrum.HiresSpectrum(targ.w[ref_order], targ.s[ref_order],
                                          targ.serr[ref_order])
    else:
        # If we already have a flattened spectrum, use the 5120-5200 A
        # region.
        targ_cut = targ.cut(5120, 5200)
        if len(targ_cut) == 0:
            # if 5120-5200 A region is missing, use the entire spectrum
            targ_cut = targ

    # Find the height of the correlation peak with each reference.
    median_peaks = []
    for i in range(len(ref_list)):
        ref = ref_list[i]
        shift_data = {}
        shift(targ_cut, ref, store=shift_data)

        # get correlation peaks
        num_sects = shift_data['order_0/num_sections']
        peaks = []
        for sect in range(num_sects):
            xcorr = shift_data['order_0/sect_{0:d}/xcorr'.format(sect)]
            if len(xcorr) > 0:
                peaks.append(max(xcorr))

        med_peak = np.median(peaks)

        median_peaks.append(med_peak)

        print("Attempting shift to spectrum {0}, ".format(ref.name) +
              "median cross-correlation peak = {0:.2f}".format(med_peak))

        if store is not None:
            store['median_peaks/{0:d}'.format(i)] = med_peak

    if store is not None:
        store['shift_reference'] = np.argmax(median_peaks)

    best_ref = ref_list[np.argmax(median_peaks)]
    print("Best reference for shifting: {0}".format(best_ref.name))

    # Now shift to the best reference
    print("Shifting entire spectrum")
    shifted = shift(targ, best_ref, store=store)

    return shifted


def shift(targ, ref, store=None, lowfilter=20):
    """Shifts the given spectrum by placing it on the same wavelength
    scale as the specified reference spectrum, then solves for shifts
    between the two spectra through cross-correlation.

    Args:
        targ (Spectrum): Target spectrum
        ref (Spectrum): Reference spectrum
        store (optional [file or dict]): h5 file or dict to record
            diagnostic data.

    Returns:
        shifted (Spectrum): Adjusted and flattened spectrum
    """
    s = np.copy(targ.s)
    serr = np.copy(targ.serr)
    w = np.copy(targ.w)
    mask = np.copy(targ.mask)
    if s.ndim == 1:
        s = np.array([s])
        serr = np.array([serr])
        w = np.array([w])
        mask = np.array([mask])

    ref = _extend_ref(ref, w[0, 0], w[-1, -1])

    # normalize each order of the target spectrum by dividing by the
    # 95th percentile
    percen_order = np.nanpercentile(s, 95, axis=1)
    s /= percen_order.reshape(-1, 1)

    # create empty 2d arrays to store each order
    s_shifted = np.asarray([[]])
    serr_shifted = np.asarray([[]])
    mask_shifted = np.asarray([[]])
    ws = np.asarray([[]])

    # create lists to store diagnostic data
    lag_data = []
    center_pix_data = []
    fit_data = []

    # length of each section in pixels
    section_length = 500

    if store is not None:
        store['num_orders'] = s.shape[0]

    # Containers for tempoerary holding of data
    s_rescaled = []
    serr_rescaled = []
    m_rescaled = []
    start_idxs = []

    # Fixed number of sections across every order
    num_sections = int(s.shape[1] / section_length) + 1

    # shift each order
    for i in range(s.shape[0]):
        ww = w[i]
        ss = s[i]
        sserr = serr[i]
        mm = mask[i]

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
        start_idxs.append(np.argmax(in_range))
        w_ref_c = ref.w[in_range]
        s_ref_c = ref.s[in_range]
        m_ref_c = ref.mask[in_range]

        # place the target spectrum on the same wavelength scale
        ss, sserr, mm = rescale_w(ss, sserr, ww, mm, w_ref_c)

        s_rescaled.append(ss)
        serr_rescaled.append(sserr)
        m_rescaled.append(mm)

        # true section length
        l_sect = int(len(ss) / num_sections)

        lags = np.empty(num_sections)
        center_pix = np.empty(num_sections)

        if store is not None:
            key = "order_{0:d}/num_sections".format(i)
            store[key] = num_sections

        for j in range(num_sections):
            # Get indices for section
            idx_min = j * l_sect
            idx_max = (j+1) * l_sect
            center_pix[j] = (j + 1/2)*l_sect

            ss_sect = ss[idx_min:idx_max]
            mm_sect = mm[idx_min:idx_max]
            s_ref_sect = s_ref_c[idx_min:idx_max]
            m_ref_sect = m_ref_c[idx_min:idx_max]

            # Don't use segments which have too many nans
            if len(ss_sect[mm_sect]) < (l_sect / 2) or \
                    len(s_ref_sect[m_ref_sect]) < (l_sect / 2):
                lag = np.nan
                lag_arr = []
                xcorr = []
            else:
                # get the shifts in pixel number
                lag, lag_arr, xcorr = solve_for_shifts(ss_sect, mm_sect,
                                                       s_ref_sect, m_ref_sect,
                                                       lowfilter=lowfilter)
            # Save results
            lags[j] = lag
            if store is not None:
                key = "order_{0:d}/sect_{1:d}/".format(i, j)
                store[key+"xcorr"] = xcorr
                store[key+"lag_arr"] = lag_arr

        # Save lag data
        lag_data.append(lags)
        center_pix_data.append(center_pix)

    lag_data = np.asarray(lag_data)
    # Compute sigma-clipped mean lags for each segment, if there are multiple
    # orders
    if s.shape[0] > 1:
        clip = 2
        for j in range(lag_data.shape[1]):
            lag_order = lag_data[:, j]
            lag_order = lag_order[~np.isnan(lag_order)]
            clipped, crit_low, crit_high = sigmaclip(lag_order, low=clip,
                                                     high=clip)

            mean_lag = np.nanmean(clipped)

            # Replace values outside the critical range and nans with mean_lag
            for i in range(lag_data.shape[0]):
                curr = lag_data[i, j]
                lag_data[i, j] = curr if curr > crit_low and curr < crit_high \
                    else mean_lag
    else:
        for j in range(lag_data.shape[1]):
            curr = lag_data[0, j]
            if j == 0:
                lag_data[0, j] = lag_data[0, j + 1] if np.isnan(curr) else curr
            elif j == (lag_data.shape[1] - 1):
                lag_data[0, j] = lag_data[0, j - 1] if np.isnan(curr) else curr
            else:
                lag_data[0, j] = np.nanmean([lag_data[0, j - 1],
                    lag_data[0, j + 1]]) if np.isnan(curr) else curr

    for i in range(s.shape[0]):
        # Restore data from previous loop
        ss = s_rescaled[i]
        sserr = serr_rescaled[i]
        mm = m_rescaled[i]
        start_idx = start_idxs[i]

        lags = lag_data[i]
        center_pix = center_pix_data[i]

        # use robust least squares to fit a line to the shifts
        # (Cauchy loss function)
        p_guess = np.array([0, 0])
        fit_res = least_squares(_linear_fit_residuals, p_guess,
                                args=(center_pix, lags), loss='cauchy')
        fit = fit_res.x
        pix_arr = np.arange(0, len(ss))
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

        # save diagnostic data
        fitted = fit[0] * center_pix + fit[1]
        fit_data.append(np.array(fitted))

    # save diagnostic data
    if store is not None:
        # convert jagged array to rectangular one
        lengths = []
        for l in lag_data:
            lengths.append(len(l))
        ml = max(lengths)

        lag_data = [utils.extend_array(l, ml) for l in lag_data]
        center_pix_data = [utils.extend_array(l, ml) for l in center_pix_data]
        fit_data = [utils.extend_array(l, ml) for l in fit_data]

        store['lag'] = np.asarray(lag_data)
        store['center_pix'] = np.asarray(center_pix_data)
        store['fit'] = np.asarray(fit_data)

    # flatten spectrum
    w_min = ws[0]
    w_max = ws[-1]
    in_range = np.asarray([True if wr > w_min and wr < w_max
                           else False for wr in ref.w])
    w_ref_trunc = ref.w[in_range]

    w_flat, s_flat, serr_flat, mask_flat = \
        flatten(ws, s_shifted, serr_shifted, mask_shifted, w_ref=w_ref_trunc)

    return spectrum.Spectrum(w_flat, s_flat, serr_flat, name=targ.name,
                             mask=mask_flat, header=targ.header,
                             attrs=targ.attrs)


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


def _extend_ref(ref, min_w, max_w):
    """Extends the reference spectrum to the given limits, assuming a constant
    delta log-lambda scale.

    Args:
        ref (spectrum.Spectrum): Reference spectrum
        min_w, max_w (float): Wavelength limts
    """
    # Delta log-lambda
    w = ref.w
    dw = np.median(np.log10(w[1:]) - np.log10(w[:-1]))

    if min_w < w[0]:
        num_points = int((np.log10(w[0]) - np.log10(min_w))/dw)
        left = np.logspace(np.log10(w[0]), np.log10(min_w), num_points,
                           base=10.0)[1:]
        # Don't forget to reverse left
        w = np.concatenate((left[::-1], w))

    if max_w > w[-1]:
        num_points = int((np.log10(max_w) - np.log10(w[-1]))/dw)
        right = np.logspace(np.log10(w[-1]), np.log10(max_w), num_points,
                            base=10.0)[1:]
        w = np.concatenate((w, right))

    if len(w) != len(ref.w):
        ref = ref.extend(w)

    return ref


def flatten(w, s, serr=None, mask=None, w_ref=None, wavlim=None):
    """Flattens a given 2-D spectrum into a 1-D array.
    Merges overlapped points by taking the mean.
    If w_ref is given, fills values that don't occur in the 2D spectrum
    with np.nan

    Args:
        w (np.ndarray): Wavelength array
        s (np.ndarray): Spectrum
        serr (np.ndarray, optional): Uncertainty in spectrum
        mask (np.ndarray, optional): Boolean mask
        w_ref (np.nadarray, optional): Reference, 1-D wavelength array
        wavlim (2-element iterable, optional): Wavelength limits

    Returns:
        w, s, serr: Wavelength, spectrum and uncertainty in spectrum
    """
    if np.shape(w) != np.shape(s):
        raise ValueError("Error: w, s not the same shape.")
    if serr is not None and np.shape(w) != np.shape(serr):
        raise ValueError("Error: w, serr not the same shape.")

    # no need to do anything if spectrum is already flat
    if np.shape(w)[0] == 1:
        w_flattened = w[0]
        s_flattened = s[0]
        serr_flattened = None if serr is None else serr[0]
        mask_flattened = None if mask is None else mask[0]
        return w_flattened, s_flattened, serr_flattened, mask_flattened

    if w_ref is None:
        w_flattened = np.unique(w)
    else:
        w_flattened = w_ref

    # create new arrays to contain spectrum
    s_flattened = np.empty_like(w_flattened)
    serr_flattened = None if serr is None else np.empty_like(w_flattened)
    mask_flattened = None if mask is None else np.empty_like(w_flattened)
    idx_max = len(w)-1
    c_idx = 0
    n_idx = 0

    for i, wl in enumerate(w_ref):
        while w[c_idx] < wl and c_idx < idx_max and not _isclose(w[c_idx], wl):
            c_idx += 1

        if c_idx >= idx_max:
            s_flattened[i] = np.nan
            if serr is not None:
                serr_flattened[i] = np.nan
            if mask is not None:
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
                    overlap = True
                    break
                elif w[n_idx] > w[c_idx] and w[n_idx] < w[n_idx-1]:
                    # lock at start of new order
                    n_idx -= 1
                    break
            if overlap:
                s_flattened[i] = (s[c_idx]+s[n_idx])/2
                if serr is not None:
                    serr_flattened[i] = (serr[c_idx]+serr[n_idx])/2
                if mask is not None:
                    mask_flattened[i] = bool(mask[c_idx]) & bool(mask[n_idx])
            else:
                s_flattened[i] = s[c_idx]
                if serr is not None:
                    serr_flattened[i] = serr[c_idx]
                if mask is not None:
                    mask_flattened[i] = mask[c_idx]
            c_idx += 1
        else:
            s_flattened[i] = np.nan
            if serr is not None:
                serr_flattened[i] = np.nan
            if mask is not None:
                mask_flattened[i] = False

    return w_flattened, s_flattened, serr_flattened, mask_flattened


def solve_for_shifts(s, mask, s_ref, mask_ref, lowfilter=20, window=1):
    """
    Solve for the pixel shifts required to align two spectra that are on the
    same wavelength scale.

    Correlates the two spectra, then fits a quadratic to the peak in order to
    solve for sub-pixel shifts.

    Args:
        s: The target spectrum array
        mask: Mask array for the target spectrum
        s_ref: The reference spectrum
        mask_ref: Mask array for the reference spectrum
        window : We fit a quadratic to pixels in [-window, window] around the peak

    Returns:
        The pixel shift, the lag and correlation data
    """
    # set masked values to nan
    s = s.copy()
    s[~mask] = np.nan
    s_ref = s_ref.copy()
    s_ref[~mask_ref] = np.nan

    # find the mean of the two spectra
    smean = np.nanmean(s)
    srefmean = np.nanmean(s_ref)
    # fill nans with mean value so they contribute nothing to correlation
    s = _fill_nans(s, smean)
    s_ref = _fill_nans(s_ref, srefmean)

    # # perform correlation
    # xcorr = np.correlate(s-smean, s_ref-srefmean, mode='full')
    # xcorr = np.nan_to_num(xcorr)
    # max_corr = np.argmax(xcorr)

    # # number of pixels
    # npix = xcorr.shape[0]
    # lag_arr = np.arange(-(npix-1)/2, (npix+1)/2, 1)

    # perform correlation
    xcorr = correlate(s-smean, s_ref-srefmean, lowfilter=lowfilter)
    max_corr = np.argmax(xcorr)
    npix = len(xcorr)
    lag_arr = np.arange(-npix/2+1, npix/2+1, 1)

    # select points around the peak and fit a quadratic
    lag_peaks = lag_arr[max_corr-window:max_corr+window+1]
    xcorr_peaks = xcorr[max_corr-window:max_corr+window+1]

    p = np.polyfit(lag_peaks, xcorr_peaks, 2)
    # peak is simply -p[1]/2p[0]
    lag = -p[1] / (2*p[0])

    return lag, lag_arr, xcorr


def correlate(a, v, lowfilter=0):
    """Custom function to perform 1-dimensional cross-correlation

    Args:
        a (np.ndarray): Input sequence
        v (np.ndarray): Input sequence
        lowfilter (int): Filter out components with wavelength above this
            number of pixels

    Returns:
        np.ndarray: Symmetric cross-correlation array
    """
    # Zero pad arrays to double length, rounded to nearest power of two
    # This is necessary to avoid effects from circular convolution
    l = max(len(a), len(v))
    # Add 1 in the exponent for doubling the length
    # Add another for integer rounding
    padded_length = 2**int(np.log2(l) + 1 + 1)
    a_padded = np.zeros(padded_length)
    a_padded[:len(a)] = a
    v_padded = np.zeros(padded_length)
    v_padded[:len(v)] = v

    # Perform fast Fourier transforms of the input sequences
    a_f = np.fft.rfft(a_padded)
    v_f = np.fft.rfft(v_padded)

    cutoff = int(padded_length/lowfilter)
    # Use a sigmoid filter
    pix = np.arange(len(a_f))
    b = expit(0.01 * (pix - cutoff))

    # Perform low frequency filtering
    a_f *= b
    v_f *= b

    # Use correlation theorem to perform a fast cross-correlation
    xcorr = np.fft.irfft(np.conj(a_f) * v_f)

    # Final correlation array is by convention reversed and reorderd
    xcorr_rearranged = np.empty(len(xcorr))
    xcorr_rearranged[:int(len(xcorr)/2)] = xcorr[int(len(xcorr)/2):]
    xcorr_rearranged[int(len(xcorr)/2):] = xcorr[:int(len(xcorr)/2)]
    xcorr_rearranged = xcorr_rearranged[::-1]

    return xcorr_rearranged


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


def shift_data_to_hdu(shift_data):
    """Saves the shift data to a BinTableHDU.

    Args:
        shift_data (dict): Shift data output from shift()
    """
    shift_data = shift_data.copy()
    col_list = []
    num_orders = shift_data.pop('num_orders')
    col_list.append(fits.Column(name='num_orders', format='J',
                                array=[num_orders]))

    num_sects = []
    for i in range(num_orders):
        num_sects.append(shift_data.pop('order_{0:d}/num_sections'
                                        .format(i)))
    col_list.append(fits.Column(name='num_sects', format='J',
                                array=num_sects))
    n_sec = max(num_sects)

    for k in ['center_pix', 'lag', 'fit']:
        col_list.append(fits.Column(name=k,
                        format='{0:d}E'.format(n_sec),
                        array=shift_data.pop(k)))

    # Save individual fit data
    for k in shift_data.keys():
        if np.ndim(shift_data[k]) == 0:
            arr = [shift_data[k]]
        else:
            arr = shift_data[k]
        col_list.append(fits.Column(name=k, format='D',
                                    array=arr))

    shift_hdu = fits.BinTableHDU.from_columns(col_list)
    shift_hdu.name = 'SHIFTDATA'

    return shift_hdu


def save_shift_to_fits(outpath, shifted, unshifted, shift_data, clobber=False):
    """Saves the complete shift data to a FITS file.

    Args:
        outpath (str): Path to save output file
        shifted (Spectrum): Shifted spectrum
        unshifted (HiresSpectrum): Raw spectrum
        shift_data (dict): Shift data
        clobber (bool): Overwrite existing file at destination
    """
    # Create primary HDU
    prihdu = fits.PrimaryHDU(header=shifted.header)
    # Save shifted spectrum
    shifted_hdu = shifted.to_hdu()

    # Save unshifted spectrum
    unshifted_hdus = unshifted.to_hdulist(primary=False)

    # Save shift data
    shift_data_hdu = shift_data_to_hdu(shift_data)

    hdulist = fits.HDUList([prihdu, shifted_hdu, shift_data_hdu] +
                           unshifted_hdus)

    hdulist.writeto(outpath, overwrite=clobber)
