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

from specmatchemp.io import pdplus
from specmatchemp import plots

class Spectrum(object):
    """Spectrum class

    Attributes:
        w (np.ndarray): Wavelength scale
        s (np.ndarray): Spectrum
        serr (None or np.ndarray): Measurement error in spectrum
        mask (np.ndarray): Boolean array with telluric line positions
        name (str): Name associted with spectrum
        header (FITS header): Header from FITS file
        attrs (dict): A dictionary of further attributes

    This object is a container for a spectrum and other
    properties which may be required.
    """
    _counter = 1
    def __init__(self, w, s=None, serr=None, mask=None, name=None, header=None, attrs={}):
        """
        Args:
            w (np.ndarray): Wavelength scale
            s (np.ndarray): Spectrum
            serr (optional [np.ndarray]): Error in spectrum
            mask (optional [np.ndarray]): Boolean array to mask out telluric lines
            name (optional [str]): Name associated with spectrum
            header (optional [FITS header]): Header from fits file
            attrs (optional [dict]): Any further attributes
        """
        self.w = w
        
        if s is None:
            self.s = np.zeros_like(w)
        else:
            self.s = s

        if serr is None:
            self.serr = np.zeros_like(serr)
        else:
            self.serr = serr

        if mask is None:
            self.mask = np.empty_like(s)
            self.mask.fill(True)
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
        return type(self)(np.copy(self.w), np.copy(self.s), np.copy(self.serr)\
            , mask=np.copy(self.mask), name=self.name, header=self.header, attrs=self.attrs.copy())


    def to_fits(self, outfile):
        """Saves the spectrum to a fits file.

        Args:
            outfile (str): Path to output file
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
            fits.Column(name='mask', format='L', array=self.mask)])

        hdulist = fits.HDUList([prihdu, tbhdu])
        hdulist.writeto(outpath, clobber=True)


    def to_hdf(self, outfile):
        """Saves the spectrum to a hdf file.

        Args:
            outfile (str or h5 file): Output path or file handle
        """
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

        outfile.create_dataset('s', data=self.s)
        outfile.create_dataset('serr', data=serr)
        outfile.create_dataset('w', data=self.w)
        outfile.create_dataset('mask', data=self.mask)

        outfile.attrs['name'] = self.name
        outfile.attrs['header'] = str(self.header)

        for key in self.attrs.keys():
            outfile.attrs[key] = self.attrs[key]

        if is_path:
            outfile.close()


    def cut(self, minw, maxw):
        """Truncate the spectrum between the given limits

        Args:
            minw (float): Minimum wavelength
            maxw (float): Maximum wavelength
        Returns:
            truncated (Spectrum): Truncated spectrum object
        """
        wavmask = (self.w >= minw) & (self.w <= maxw)

        w_trunc = self.w[wavmask]
        s_trunc = self.s[wavmask]
        serr_trunc = None if self.serr is None else self.serr[wavmask]
        mask_trunc = None if self.mask is None else self.mask[wavmask]

        return type(self)(w_trunc, s_trunc, serr_trunc, mask_trunc, self.name, self.attrs)


    def plot(self, offset=0, label='_nolegend_', plt_kw={'color':'RoyalBlue'}, showmask=False,\
            text='', text_kw={}):
        """Plots the spectrum

        Args:
            offset (optional [float]): Vertical offset of the spectrum
            label (optional [str]): Label of spectrum (appears in plt.legend)
            plt_kw (optional [dict]): Keyword arguments to pass to plt.plot
            text (optional [str]): String to label the spectrum
            text_kw (optional [dict]): Keyword arguments to pass to plt.text
        """
        plt.plot(self.w, self.s+offset, '-', label=label, **plt_kw)
        if len(text) > 0:
            plots.annotate_spectrum(text, spec_offset=offset, text_kw=text_kw)

        if showmask:
            ax = plt.gca()
            # get list of masked regions
            regions = self._convert_mask_to_regions()
            for reg in regions:
                ax.add_patch(patches.Rectangle((reg[0],offset), reg[1]-reg[0], offset+1,\
                    ec='none', fc='gray', alpha=0.3))

        plt.grid(True)
        plt.xlabel('Wavelength (Angstroms)')
        plt.ylabel('Normalized Flux (Arbitrary Offset)')


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
                    l.append((w[start],w[end]))
            if ismasked:
                end = len(self.mask)
                l.append((w[start],w[end]))
            return l
        ###

        l = []
        if self.w.ndim == 1:
            l = get_regions(self.mask, self.w)
        else:
            for i in range(self.w.shape[0]):
                l.append(get_regions(self.mask[i], self.w[i]))

        return l


class HiresSpectrum(Spectrum):
    """Spectrum class for raw HIRES spectra.

    Is read by the read_hires
        
    Attributes:
        w (np.ndarray): Wavelength scale
        s (np.ndarray): Spectrum
        serr (None or np.ndarray): Measurement error in spectrum
        mask (np.ndarray): Boolean array with telluric line positions
        name (str): Name associted with spectrum
        header (FITS header): Header from FITS file
        attrs (dict): A dictionary of further attributes
    """
    def __init__(self, w, s, serr=None, mask=None, mask_table=None, name=None, header=None, attrs={}):
        """
        Args:
            w (np.ndarray): Wavelength scale
            s (np.ndarray): Spectrum
            serr (optional [np.ndarray]): Error in spectrum
            mask (optional [np.ndarray]): Boolean array to mask out telluric lines
            mask_table (optional [pd.DataFrame]): Table containing masked regions
            name (optional [str]): Name associated with spectrum
            header (optional [FITS header]): Header from fits file
            attrs (optional [dict]): Any further attributes
        """
        self.mask_table = mask_table
        if mask_table is not None and mask is None:
            mask = np.empty_like(s)
            mask.fill(True)
            for order in range(len(mask)):
                mask_table_cut = mask_table.query('order == {0:d}'.format(order))
                for n, row in mask_table_cut.iterrows():
                    start = row['minpix']
                    end = row['maxpix']
                    mask[order,start:end] = False

        super(HiresSpectrum, self).__init__(w, s, serr, mask=mask, name=name, header=header, attrs=attrs)

    def plot(self, offset=0, normalize=False, label='_nolegend_', plt_kw={}, showmask=False):
        """Plots the spectrum

        Args:
            offset (optional [float]): Vertical offset of the spectrum
            label (optional [str]): Label of spectrum (appears in plt.legend)
            plt_kw (optional [dict]): Keyword arguments to pass to plt.plot
        """
        if normalize:
            plt.plot(self.w.T, self.s.T/np.percentile(self.s, 95, axis=1)+offset, '-', label=label, **plt_kw)
        else:
            plt.plot(self.w.T, self.s.T+offset, '-', label=label, **plt_kw)

        if showmask:
            ax = plt.gca()
            ylim = ax.get_ylim()
            # get list of masked regions
            regions = self._convert_mask_to_regions()
            for order in range(len(self.w)):
                for reg in regions[order]:
                    ax.add_patch(patches.Rectangle((reg[0],ylim[0]), reg[1]-reg[0], ylim[1]-ylim[0],\
                        ec='none', fc='gray', alpha=0.3))

        plt.grid(True)
        plt.xlabel('Wavelength (Angstroms)')
        plt.ylabel('Flux')


def read_fits(infile, wavlim=None):
    """Reads a spectrum from a fits file

    Args:
        infile (str): Path to input fits file
        wavlim (optional [tuple]): Wavelength limits to read
    Returns:
        spec (Spectrum): Spectrum object
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
    header = hdu[0].header

    hdu.close()

    if wavlim is None:
        return Spectrum(w, s, serr, mask=mask, header=header)
    else:
        return Spectrum(w, s, serr, mask=mask, header=header).cut(*wavlim)


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

    mask_table = None
    if maskfile is not None:
        mask_table = pd.read_csv(maskfile)
        # Get HIRES chip
        chip = os.path.basename(infile)[0:2]
        mask_table = mask_table[mask_table.chip.str.contains(chip)]


    return HiresSpectrum(w, s, serr, mask_table=mask_table, header=header)

def read_hdf(infile):
    """Reads a spectrum from a hdf file

    Args:
        infile (str or h5 file): Input path or file handle
    Returns:
        spec (Spectrum): Spectrum object
    """
    is_path = False
    if isinstance(infile, str):
        infile = h5py.File(infile, 'r')
        is_path = True

    s = infile['s'][:]
    serr = infile['serr'][:]
    w = infile['w'][:]
    if 'mask' in infile.keys():
        mask = infile['mask'][:]
    else:
        mask = None
    
    attrs = dict(infile.attrs)
    if 'name' in attrs.keys():
        name = attrs.pop('name')
    else:
        name=None

    if 'header' in attrs.keys():
        header = attrs.pop('header')
    else:
        header=None

    if is_path:
        infile.close()

    return Spectrum(w, s, serr, mask=mask, name=name, header=header, attrs=attrs)

