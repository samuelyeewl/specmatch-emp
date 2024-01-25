"""
@filename spectrum.py

Defines a spectrum object for use in SpecMatch-Emp
"""

import os
from warnings import catch_warnings, filterwarnings
import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy.io import fits
from scipy.signal.windows import gaussian
from scipy.ndimage import convolve1d

from specmatchemp import plots


class Spectrum(object):
    """Spectrum class.

    This object is a container for a spectrum and its properties.

    Args:
        w (np.ndarray): Wavelength scale
        s (np.ndarray): Spectrum
        serr (np.ndarray, optional): Error in spectrum
        mask (np.ndarray, optional): Boolean array to mask out telluric lines
        name (str, optional): Name associated with spectrum
        header (FITS header, optional): Header from fits file
        attrs (dict, optional): Any further attributes

    Attributes:
        w (np.ndarray): Wavelength scale
        s (np.ndarray): Spectrum
        serr (np.ndarray): Measurement error in spectrum
        mask (np.ndarray): Boolean array with telluric line positions
        name (str): Name associted with spectrum
        header (FITS header): Header from FITS file
        attrs (dict): A dictionary of further attributes
    """
    _counter = 1

    def __init__(self, w, s, serr=None, mask=None, name=None, header=None,
                 attrs={}):
        # Wavelength and spectrum are both required
        self.w = w
        self.s = s

        if serr is None:
            self.serr = np.zeros_like(w)
        else:
            self.serr = serr

        if mask is None:
            self.mask = np.empty_like(s).astype(bool)
            self.mask.fill(True)
        elif isinstance(mask, np.ndarray):
            self.mask = mask.astype(bool)
        else:
            self.mask = mask

        if name is None:
            self.name = "Spectrum {0:d}".format(Spectrum._counter)
            Spectrum._counter += 1
        else:
            self.name = name

        self.attrs = attrs
        self.header = header

    def copy(self):
        """Returns a deep copy of the Spectrum object
        """
        return type(self)(np.copy(self.w), np.copy(self.s), np.copy(self.serr),
                          mask=np.copy(self.mask), name=self.name,
                          header=self.header, attrs=self.attrs.copy())

    def to_fits(self, outpath, clobber=True):
        """Saves the spectrum to a fits file.

        Args:
            outpath (str): Path to output file
        """
        if self.header is None:
            self.header = fits.Header()
        prihdu = fits.PrimaryHDU(header=self.header)

        # If no error is provided
        if self.serr is None:
            serr = np.empty_like(self.s)
            serr.fill(np.nan)
        else:
            serr = self.serr

        tbhdu = fits.BinTableHDU.from_columns(
            [fits.Column(name='s', format='D', array=self.s),
             fits.Column(name='w', format='D', array=self.w),
             fits.Column(name='serr', format='D', array=serr),
             fits.Column(name='mask', format='B', array=self.mask)])

        hdulist = fits.HDUList([prihdu, tbhdu])
        hdulist.writeto(outpath, overwrite=clobber)

    def to_hdu(self):
        """Creates a fits.BinTableHDU object from spectrum data

        Returns:
            fits.BinTableHDU:
                A binary table HDU object
        """
        hdu = fits.BinTableHDU.from_columns(
            [fits.Column(name='s', format='D', array=self.s),
             fits.Column(name='w', format='D', array=self.w),
             fits.Column(name='serr', format='D', array=self.serr),
             fits.Column(name='mask', format='L', array=self.mask)])

        # add metadata
        hdu.header['NAME'] = self.name
        for k in self.attrs.keys():
            hdu.header[k[0:8]] = self.attrs[k]

        return hdu

    def to_hdf(self, outfile, suffix=""):
        """Saves the spectrum to a hdf file.

        Args:
            outfile (str or h5 file): Output path or file handle
            suffix (str, optional): Suffix to append to h5 field names
        """
        # Allow either a string or h5 file object ot be passed.
        is_path = False
        if isinstance(outfile, str):
            outfile = h5py.File(outfile, 'w')
            is_path = True

        # If no error is provided
        if self.serr is None:
            serr = np.empty_like(self.s)
            serr.fill(np.nan)
        else:
            serr = self.serr

        outfile.create_dataset('s' + suffix, data=self.s)
        outfile.create_dataset('serr' + suffix, data=serr)
        outfile.create_dataset('w' + suffix, data=self.w)
        outfile.create_dataset('mask' + suffix, data=self.mask)

        outfile.attrs['name' + suffix] = self.name
        outfile.attrs['header' + suffix] = str(self.header)

        for key in self.attrs.keys():
            outfile.attrs[key + suffix] = self.attrs[key]

        if is_path:
            outfile.close()

    def cut(self, minw, maxw):
        """Truncate the spectrum between the given limits

        Args:
            minw (float): Minimum wavelength
            maxw (float): Maximum wavelength
        Returns:
            Spectrum: Truncated spectrum object
        """
        wavmask = (self.w >= minw) & (self.w <= maxw)

        w_trunc = self.w[wavmask]
        s_trunc = self.s[wavmask]

        serr_trunc = None if self.serr is None else self.serr[wavmask]
        mask_trunc = None if self.mask is None else self.mask[wavmask]

        return type(self)(w_trunc, s_trunc, serr_trunc, mask_trunc,
                          name=self.name, attrs=self.attrs.copy())

    def extend(self, w):
        """Extend the spectrum onto a new wavelength scale, placing np.nan
        for points which do not exist on the old scale.

        Args:
            w (np.ndarray): New wavlength array.
        """
        wavmap = np.searchsorted(w, self.w)
        s_extend = np.empty_like(w, dtype=np.float)
        s_extend.fill(np.nan)
        s_extend[wavmap] = self.s
        serr_extend = np.empty_like(w, dtype=np.float)
        serr_extend.fill(np.nan)
        serr_extend[wavmap] = self.serr
        mask_extend = np.empty_like(w, dtype=bool)
        mask_extend.fill(False)
        mask_extend[wavmap] = self.mask

        return type(self)(w, s_extend, serr_extend, mask_extend,
                          name=self.name, attrs=self.attrs.copy())

    def rescale(self, w):
        """Puts the spectrum onto a new wavelength scale, interpolating between
        points if necessary. Places np.nan for points beyond existing range.

        Args:
            w (np.ndarray): New wavelength scale.
        """
        if self.on_scale(w):
            return self

        snew = np.interp(w, self.w, self.s, left=np.nan, right=np.nan)
        serrnew = np.interp(w, self.w, self.serr, left=np.nan, right=np.nan)
        masknew = np.interp(w, self.w, self.mask, left=0, right=0).astype(bool)
        return type(self)(w, snew, serrnew, masknew,
                          name=self.name, attrs=self.attrs.copy())

    def on_scale(self, w):
        """Checks if spectrum is on the given wavelength scale

        Args:
            w (np.ndarray): Wavelength scale to check against
        """
        return np.allclose(w, self.w, rtol=1e-8)

    def wavlim(self):
        """Gets the wavelength range of the spectrum"""
        return (self.w[0], self.w[-1])

    def plot(self, wavlim='all', offset=0, label='_nolegend_', showmask=False,
             normalize=True, plt_kw={'color': 'RoyalBlue'}, text='',
             text_kw={}):
        """Plots the spectrum.

        Args:
            offset (float, optional): Vertical offset of the spectrum
            label (str, optional): Label of spectrum (appears in plt.legend)
            showmask (bool, optional): Whether to highlight the telluric mask
                in the plot. Defaults to False.
            plt_kw (dict, optional): Keyword arguments to pass to plt.plot
            text (str, optional): String to label the spectrum
            text_kw (dict, optional): Keyword arguments to pass to plt.text
        """
        if normalize:
            plt.plot(self.w, self.s / np.nanpercentile(self.s, 95) + offset,
                     '-', label=label, **plt_kw)
        else:
            plt.plot(self.w, self.s + offset, '-', label=label, **plt_kw)
        if len(text) > 0:
            plots.annotate_spectrum(text, spec_offset=offset, text_kw=text_kw)

        if showmask:
            ax = plt.gca()
            # get list of masked regions
            regions = self._convert_mask_to_regions()
            for reg in regions:
                ax.add_patch(patches.Rectangle((reg[0], offset),
                                               reg[1] - reg[0], offset + 1,
                                               ec='none', fc='gray',
                                               alpha=0.3))

        plt.grid(True)
        if wavlim != 'all' and isinstance(wavlim, tuple):
            plt.xlim(wavlim)
        plt.xlabel('Wavelength (Angstroms)')
        plt.ylabel('Normalized Flux (Arbitrary Offset)')

    def snr(self):
        return np.nanpercentile(1/self.serr, 90)

    def _convert_mask_to_regions(self):
        """Converts a boolean mask into a list of regions
        """
        # p is FALSE for values to be masked out
        # inner function to handle flat and stacked spectra
        def get_regions(mask, w):
            ismasked = False
            start = 0
            end = 0
            l = []
            for i, p in enumerate(mask):
                if not ismasked and not p:
                    ismasked = True
                    start = i
                elif ismasked and p:
                    ismasked = False
                    end = i
                    l.append((w[start], w[end]))
            if ismasked:
                end = len(self.mask)
                l.append((w[start], w[end]))
            return l
        ###

        l = []
        if self.w.ndim == 1:
            l = get_regions(self.mask, self.w)
        else:
            for i in range(self.w.shape[0]):
                l.append(get_regions(self.mask[i], self.w[i]))

        return l

    def __len__(self):
        """Get number of elements in spectrum.
        """
        return self.w.size

    def broaden_specres(self, spec_res, kernel_npoints=151,
                        convolve_mode='nearest', **convolve_kwargs):
        """Broaden spectrum to a target spectral resolution.
        """
        mean_w = np.mean(self.w)
        gaussian_fwhm = mean_w / spec_res
        gaussian_sigma = gaussian_fwhm / 2.355
        pixels_per_A = len(self.w) / (self.w[-1] - self.w[0])
        gaussian_sigma_pix = gaussian_sigma * pixels_per_A

        kernel = gaussian(kernel_npoints, gaussian_sigma_pix)

        broadened_spec = self.copy()
        broadened_spec.s = convolve1d(self.s, kernel,
                                      mode=convolve_mode, **convolve_kwargs)

        kernel_sq = kernel **2
        kernel_sq /= np.sum(kernel_sq)
        broadened_spec.serr = np.sqrt(convolve1d(self.serr**2, kernel_sq,
                                                 mode=convolve_mode, **convolve_kwargs))

        return broadened_spec

    @staticmethod
    def combine_spectra(spectra, w_ref, name=None, prefixes=None):
        """Combines several spectra into one spectrum object.

        Combine multiple spectrum objects into a single spectrum.
        Places np.nan for the gaps between the spectra, and averages
        overlapping portions.

        Args:
            spectra (iterable): Iterable of spectrum objects
            w_ref (np.ndarray): Reference wavelength scale. This scale is used
                to fill in the gaps between the spectra.
            name (str, optional): Name for the flattend spectrum
            prefixes (str, optional): Prefixes for attributes
        """
        for spec in spectra:
            if type(spec) is not Spectrum:
                raise TypeError("combine_spectra can only be used on Spectrum")

        global_min = float('inf')
        global_max = 0.
        rescaled_spectra = []
        for spec in spectra:
            # Check wavelength scales
            min_w = spec.w[0]
            max_w = spec.w[-1]
            global_min = min_w if min_w < global_min else global_min
            global_max = max_w if max_w > global_max else global_max

            w_trunc = w_ref[(w_ref >= min_w) & (w_ref <= max_w)]
            spec = spec.rescale(w_trunc)
            rescaled_spectra.append(spec)

        # Get global wavelength array
        w = w_ref[(w_ref >= global_min) & (w_ref <= global_max)]
        s_stacked = []
        serr_stacked = []
        mask_stacked = []
        for spec in rescaled_spectra:
            # Extend spectra
            spec = spec.extend(w)
            s_stacked.append(spec.s)
            serr_stacked.append(spec.serr)
            mask_stacked.append(spec.mask)

        s_stacked = np.asarray(s_stacked)
        serr_stacked = np.asarray(serr_stacked)
        mask_stacked = np.asarray(mask_stacked)

        # Flatten spectra
        with catch_warnings():
            filterwarnings("ignore", message="Mean of empty slice")
            # Weighted averages
            weights = 1/(serr_stacked**2)
            inv_var = np.nansum(weights, axis=0)
            s = np.nansum(s_stacked * weights, axis=0) / inv_var
            serr = np.sqrt(inv_var)
            serr[~np.isfinite(serr)] = np.nan
            mask = np.any(mask_stacked, axis=0)

        # Create spectrum
        if name is None:
            name = spectra[0].name
        if prefixes is None:
            prefixes = [str(n) for n in np.arange(len(spectra))]
        attrs = {}
        for (i, spec) in enumerate(spectra):
            for k, v in spec.attrs.items():
                attrs[prefixes[i] + '_' + k] = v

        return Spectrum(w, s, serr, mask, name=name, header=None, attrs=attrs)


class EchelleSpectrum(Spectrum):
    """2D Echelle Spectrum

    Args:
        w (np.ndarray): Wavelength scale
        s (np.ndarray): Spectrum
        serr (np.ndarray, optional): Error in spectrum
        mask (np.ndarray, optional): Boolean array to mask out telluric lines
        name (str, optional): Name associated with spectrum
        header (FITS header, optional): Header from fits file
        attrs (dict, optional): Any further attributes

    Attributes:
        w (np.ndarray): Wavelength scale
        s (np.ndarray): Spectrum
        serr (np.ndarray): Measurement error in spectrum
        mask (np.ndarray): Boolean array with telluric line positions
        name (str): Name associted with spectrum
        header (FITS header): Header from FITS file
        attrs (dict): A dictionary of further attributes
    """
    def __getitem__(self, idx):
        """
        Get the given order in the spectrum.
        """
        order_w = self.w[idx]
        order_s = self.s[idx]
        order_serr = self.serr[idx]
        order_mask = self.mask[idx]
        return Spectrum(order_w, order_s, order_serr, order_mask,
                        name=self.name + f' Order {idx}', header=self.header,
                        attrs=self.attrs.copy())

    def cut(self, minw, maxw):
        """Truncate the spectrum between the given limits, returning the order
        that maximizes the length.

        Returns:
            Spectrum: Truncated spectrum object.
        """
        order_wavmask = [
            ((order_w >= minw) & (order_w <= maxw))
            for order_w in self.w
        ]
        order_wavcount = [m.sum() for m in order_wavmask]
        max_order = np.argmax(order_wavcount)

        order_w = self.w[max_order]
        order_s = self.s[max_order]
        order_serr = self.serr[max_order]
        order_mask = self.mask[max_order]
        wavmask = order_wavmask[max_order]

        w_trunc = order_w[wavmask]
        s_trunc = order_s[wavmask]
        serr_trunc = None if order_serr is None else order_serr[wavmask]
        mask_trunc = None if order_mask is None else order_mask[wavmask]

        return Spectrum(w_trunc, s_trunc, serr_trunc, mask_trunc,
                        name=self.name + f' Order {max_order}', header=self.header,
                        attrs=self.attrs.copy())

    @classmethod
    def from_speclist(cls, spectra, name=None, header=None, attrs=None):
        """Create EchelleSpectrum class from a list of spectra
        """
        w = [spec.w for spec in spectra]
        s = [spec.s for spec in spectra]
        serr = [spec.serr for spec in spectra]
        mask = [spec.mask for spec in spectra]

        if name is None:
            name = spectra[0].name
        if header is None:
            header = spectra[0].header
        if attrs is None:
            attrs = spectra[0].attrs

        return cls(w, s, serr, mask, name=name, header=header, attrs=attrs)


class HiresSpectrum(Spectrum):
    """Spectrum class for raw HIRES spectra.

    Is read by the read_hires_fits function.

    Attributes:
        w (np.ndarray): Wavelength scale
        s (np.ndarray): Spectrum
        serr (None or np.ndarray): Measurement error in spectrum
        mask (np.ndarray): Boolean array with telluric line positions
        name (str): Name associted with spectrum
        header (FITS header): Header from FITS file
        attrs (dict): A dictionary of further attributes

    Args:
        w (np.ndarray): Wavelength scale
        s (np.ndarray): Spectrum
        serr (np.ndarray, optional): Error in spectrum
        mask (np.ndarray, optional): Boolean array to mask out
            telluric lines
        mask_table (pd.DataFrame, optional): Table containing
            masked regions
        name (str, optional): Name associated with spectrum
        header (FITS header, optional): Header from fits file
        attrs (dict, optional): Any further attributes
    """
    def __init__(self, w, s, serr=None, mask=None, mask_table=None, name=None,
                 header=None, attrs={}):
        self.mask_table = mask_table
        if mask_table is not None and mask is None:
            mask = np.empty_like(s).astype(bool)
            mask.fill(True)
            for order in range(len(mask)):
                mask_table_cut = mask_table.query('order == {0:d}'
                                                  .format(order))
                for n, row in mask_table_cut.iterrows():
                    start = row['minpix']
                    end = row['maxpix']
                    mask[order, start:end] = False

        super(HiresSpectrum, self).__init__(w, s, serr, mask=mask, name=name,
                                            header=header, attrs=attrs)

    def plot(self, wavlim='all', offset=0, label='_nolegend_', normalize=False,
             showmask=False, plt_kw={'color': 'RoyalBlue'},
             text='', text_kw={}):
        """Plots the spectrum

        Args:
            wavlim (optional [tuple]): Wavelength limit to plot
            offset (optional [float]): Vertical offset of the spectrum
            label (optional [str]): Label of spectrum (appears in plt.legend)
            normalize (optional [bool]): Whether to normalize the spectrum.
                Defaults to False
            showmask (bool, optional): Whether to highlight the telluric mask
                in the plot. Defaults to False.
            plt_kw (optional [dict]): Keyword arguments to pass to plt.plot
            text (str, optional): String to label the spectrum
            text_kw (dict, optional): Keyword arguments to pass to plt.text
        """
        if self.w.ndim > 1:
            if normalize:
                plt.plot(self.w.T,
                         self.s.T / np.nanpercentile(self.s, 95, axis=1) + offset,
                         '-', label=label, **plt_kw)
            else:
                plt.plot(self.w.T, self.s.T + offset, '-', label=label,
                         **plt_kw)

            if showmask:
                ax = plt.gca()
                ylim = ax.get_ylim()
                # get list of masked regions
                regions = self._convert_mask_to_regions()
                for order in range(len(self.w)):
                    for reg in regions[order]:
                        ax.add_patch(patches.Rectangle((reg[0], ylim[0]),
                            reg[1] - reg[0], ylim[1] - ylim[0],
                            ec='none', fc='gray', alpha=0.3))
        else:
            if normalize:
                plt.plot(self.w, self.s / np.nanpercentile(self.s, 95) + offset,
                         '-', label=label, **plt_kw)
            else:
                plt.plot(self.w, self.s + offset, '-', label=label, **plt_kw)

        if len(text) > 0:
            plots.annotate_spectrum(text, spec_offset=offset, text_kw=text_kw)

        plt.grid(True)
        if wavlim != 'all' and isinstance(wavlim, tuple):
            plt.xlim(wavlim)
        plt.xlabel('Wavelength (Angstroms)')
        plt.ylabel('Flux')

    def to_hdulist(self, primary=True):
        """Creates a list of fits.ImageHDU to store HIRES spectrum.

        Args:
            primary (bool): True if first HDU should be a primary HDU
        Returns:
            list of fits.ImageHDU:
                l[0] = spectrum
                l[1] = spectrum error
                l[2] = wavelength
        """
        if primary:
            s_hdu = fits.PrimaryHDU(data=self.s, header=self.header)
        else:
            s_hdu = fits.ImageHDU(data=self.s, header=self.header)
        serr_hdu = fits.ImageHDU(data=self.serr)
        w_hdu = fits.ImageHDU(data=self.w)

        return fits.HDUList([s_hdu, serr_hdu, w_hdu])

    def to_hires_fits(self, outfile, clobber=False):
        hdulist = self.to_hdulist()
        hdulist.writeto(outfile, overwrite=clobber)

    @staticmethod
    def combine_spectra(spectra, name=None, prefixes=None):
        """Combines several raw HIRES spectra into one HIRES spectrum object.

        Args:
            spectra (iterable): Iterable of spectrum objects
            name (str, optional): Name for the flattend spectrum
            prefixes (str, optional): Prefixes for attributes
        """
        for spec in spectra:
            if type(spec) is not HiresSpectrum:
                raise TypeError("combine_spectra can only be used on " +
                                "HiresSpectrum")

        w_stacked = []
        s_stacked = []
        serr_stacked = []
        mask_stacked = []
        # Simply stack spectra
        for spec in spectra:
            w_stacked.append(spec.w)
            s_stacked.append(spec.s)
            serr_stacked.append(spec.serr)
            mask_stacked.append(spec.mask)

        w_stacked = np.concatenate(w_stacked)
        s_stacked = np.concatenate(s_stacked)
        serr_stacked = np.concatenate(serr_stacked)
        mask_stacked = np.concatenate(mask_stacked)

        # Create spectrum
        if name is None:
            name = spectra[0].name
        if prefixes is None:
            prefixes = [str(n) for n in np.arange(len(spectra))]
        attrs = {}
        for (i, spec) in enumerate(spectra):
            for k, v in spec.attrs.items():
                attrs[prefixes[i] + '_' + k] = v

        return HiresSpectrum(w_stacked, s_stacked, serr_stacked, mask_stacked,
                             name=name, header=None, attrs=attrs)


def read_fits(infile, wavlim=None):
    """Reads a spectrum from a fits file

    Args:
        infile (str): Path to input fits file
        wavlim (tuple, optional): Wavelength limits to read
    Returns:
        Spectrum: Spectrum object
    """
    # Exception will be raised by fits.open if file does not exist
    hdu = fits.open(infile)

    data = hdu[1].data
    s = data['s']
    serr = data['serr']
    w = data['w']
    if 'mask' in data.dtype.names:
        mask = data['mask']
    else:
        mask = None
    name = os.path.splitext(os.path.basename(infile))[0]
    header = hdu[0].header

    hdu.close()

    if wavlim is None:
        return Spectrum(w, s, serr, mask=mask, name=name, header=header)
    else:
        return Spectrum(w, s, serr, mask=mask, name=name, header=header)\
            .cut(*wavlim)


def read_hires_fits(infile, maskfile=None):
    """Reads a spectrum from a fits file that was produced by HIRES.

    Args:
        infile (str): Path to input fits file
        maskfile (optional [str]): Path to file containing telluric
            lines mask.
    Returns:
        spec (HiresSpectrum): Spectrum object
    """
    hdu = fits.open(infile)
    s = hdu[0].data
    serr = hdu[1].data
    w = hdu[2].data
    header = hdu[0].header
    hdu.close()

    # Convert zeros to nans
    s[np.isclose(s, 0.0)] = np.nan

    name = os.path.splitext(os.path.basename(infile))[0]

    mask_table = None
    if maskfile is not None:
        mask_table = pd.read_csv(maskfile)
        # Get HIRES chip
        chip = os.path.basename(infile)[0:2]
        mask_table = mask_table[mask_table.chip.str.contains(chip)]

    return HiresSpectrum(w, s, serr, name=name, mask_table=mask_table,
                         header=header)


def read_apf_fits(infile, wavfile, maskfile=None):
    """Reads a spectrum from a fits file that was produced by HIRES.

    Args:
        infile (str): Path to input fits file
        maskfile (optional [str]): Path to file containing telluric
            lines mask.
    Returns:
        spec (HiresSpectrum): Spectrum object
    """
    hdu = fits.open(infile)
    s = hdu[0].data
    serr = np.sqrt(s)
    w = fits.getdata(wavfile)
    header = hdu[0].header
    hdu.close()

    # Convert zeros to nans
    s[np.isclose(s, 0.0)] = np.nan

    name = os.path.splitext(os.path.basename(infile))[0]

    mask_table = None
    if maskfile is not None:
        mask_table = pd.read_csv(maskfile)
        # Get HIRES chip
        chip = os.path.basename(infile)[0:2]
        mask_table = mask_table[mask_table.chip.str.contains(chip)]

    return HiresSpectrum(w, s, serr, name=name, mask_table=mask_table,
                         header=header)



def read_hdf(infile, suffix=""):
    """Reads a spectrum from a hdf file

    Args:
        infile (str or h5 file): Input path or file handle
        suffix (str, optional): Suffix on h5 keys
    Returns:
        spec (Spectrum): Spectrum object
    """
    is_path = False
    if isinstance(infile, str):
        infile = h5py.File(infile, 'r')
        is_path = True

    s = infile['s' + suffix][:]
    serr = infile['serr' + suffix][:]
    w = infile['w' + suffix][:]
    if ('mask' + suffix) in infile.keys():
        mask = infile['mask' + suffix][:]
    else:
        mask = None

    attrs = dict(infile.attrs)
    if ('name' + suffix) in attrs.keys():
        name = attrs.pop('name' + suffix)
    else:
        name = None

    if ('header' + suffix) in attrs.keys():
        header = attrs.pop('header' + suffix)
    else:
        header = None

    if len(suffix) > 0:
        # strip suffixes from remaining keys
        for k in attrs.keys():
            attr = attrs.pop(k)
            attrs[k[:-len(suffix)]] = attr

    if is_path:
        infile.close()

    if w.ndim > 1:
        return HiresSpectrum(w, s, serr, mask=mask, name=name,
                             header=header, attrs=attrs)
    else:
        return Spectrum(w, s, serr, mask=mask, name=name,
                        header=header, attrs=attrs)
