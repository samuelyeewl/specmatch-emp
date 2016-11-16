"""
@filename spectrum.py

Defines a spectrum object for use in SpecMatch-Emp
"""

import os
import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy.io import fits

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
        else:
            self.mask = mask.astype(bool)

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
        hdulist.writeto(outpath, clobber=clobber)

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
        return np.allclose(w, self.w)

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
                         self.s.T / np.percentile(self.s, 95, axis=1) + offset,
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
                plt.plot(self.w, self.s / np.percentile(self.s, 95) + offset,
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
        hdulist.writeto(outfile, clobber=clobber)


def read_fits(infile, wavlim=None):
    """Reads a spectrum from a fits file

    Args:
        infile (str): Path to input fits file
        wavlim (tuple, optional): Wavelength limits to read
    Returns:
        Spectrum: Spectrum object
    """
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
