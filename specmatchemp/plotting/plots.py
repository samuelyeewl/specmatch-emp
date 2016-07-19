"""
@filename plotting.py

Helper functions to plot various data from SpecMatch-Emp
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import os
import lmfit

from specmatchemp import library
from specmatchemp import match
from specmatchemp.io import specmatchio

UNSHIFTED_PATH = '/Users/samuel/Dropbox/SpecMatch-Emp/spectra/iodfitsdb/{0}.fits'

############################## Helper functions ################################
def reverse_x():
    plt.xlim(plt.xlim()[::-1])

def reverse_y():
    plt.ylim(plt.ylim()[::-1])

######################### Library and Spectrum plots ###########################
def plot_library_params(lib, param_x, param_y, grouped=False, plt_kw={}):
    """Plot a H-R diagram from the library

    Args:
        lib (library.Library): The library object
            or (pd.DataFrame): The library params dataframe
        param_x (str): Parameter to be plot on the x-axis
        param_y (str): Parameter to be plot on the y-axis
        grouped (bool): (optional): Whether to group library by source catalog
    """
    if type(lib) is library.Library:
        params = lib.library_params
    elif type(lib) is pd.DataFrame:
        params = lib
    else:
        raise TypeError

    x = param_x
    y = param_y
    assert x in params.columns, "{0} not in library_params".format(x)
    assert y in params.columns, "{0} not in library_params".format(y)
    
    if grouped:
        g = params.groupby('source')
        for source in g.groups:
            cut = params.ix[g.groups[source]]
            plt.plot(cut[x], cut[y], '.', label=source)
    else:
        plt.plot(params[x], params[y], '.', **plt_kw)

def plot_hires_spectrum(filename, wavlim=None, label=None, offset=0):
    """Plot a HIRES spectrum within a given wavelength range

    Args:
        filename (str): FITS file containing HIRES spectrum
        wavlim (2-element iterable): (optional) Wavelength range
        label (str): (optional) Label for this specturm
        offset (float): (optional) Vertical offset of spectrum
    """
    w, s, serr, hdr = specmatchio.read_hires_spectrum(filename)
    w = w.reshape(-1)
    s = s.reshape(-1)
    if wavlim is not None:
        w, s, serr = specmatchio.truncate_spectrum(wavlim, w, s)
    percen = np.percentile(s,95)
    s /= percen
    plt.plot(w, s+offset, label=label)

def plot_standard_spectrum(filename, wavlim=None, label=None, offset=0):
    """Plot a standard spectrum saved by specmatchemp.io.io within
    a given wavelength range

    Args:
        filename (str): FITS file containing spectrum
        wavlim (2-element iterable): (optional) Wavelength range
        label (str): (optional) Label for this specturm
        offset (float): (optional) Vertical offset of spectrum
    """
    w, s, serr, hdr = specmatchio.read_standard_spectrum(filename)
    if wavlim is not None:
        w, s, serr = specmatchio.truncate_spectrum(wavlim, w, s)
    percen = np.percentile(s,95)
    s /= percen
    plt.plot(w, s+offset, label=label)

def plot_library_spectrum(lib, lib_index, wavlim=None, offset=0, plt_kw={}):
    """Plot a spectrum from the library.

    Args:
        lib (library.Library): The library object
        lib_index (int): Index of spectrum in library
        wavlim (2-element iterable): (optional) Wavelength range
        label (str): (optional) Label for this specturm
        offset (float): (optional) Vertical offset of spectrum
    """
    if lib_index not in lib:
        print("Index {0} not in library.".format(lib_index))
        return

    s = lib.library_spectra[lib_index, 0]
    serr = lib.library_spectra[lib_index,1]
    w = lib.wav
    if wavlim is not None:
        w, s, serr = specmatchio.truncate_spectrum(wavlim, w, s, serr)
    plt.plot(w, s+offset, **plt_kw)
######################### Library and Spectrum plots ###########################


################################# Shift plots ##################################
def shifted_spectrum_plot(lib, lib_index, ref, wavlim=None, offset=True):
    """Plot the shifted and unshifted spectra against the reference

    Args:
        lib (library.Library): The library object
        lib_index (int): Index of spectrum in library
        ref: Either path to a FITS file containing a standard spectrum, or
            a library index to the spectrum used as a reference.
        wavlim (2-element iterable): (optional) Wavelength range
        offset (bool): Whether to plot the spectra offset from each other.
    """
    # find reference spectrum
    if os.path.isfile(ref):
        w_ref, s_ref, serr_ref, hdr_ref = specmatchio.read_standard_spectrum(ref, wavlim)
    else:
        pattern = '^'+ref+'$'
        row = lib.library_params[lib.library_params.lib_obs.str.match(pattern)]
        if row.empty:
            print("Could not find observation {0} in library".format(ref))
            return
        ref_idx = row.iloc[0].lib_index
        s_ref = lib.library_spectra[ref_idx,0]
        serr_ref = lib.library_spectra[ref_idx,1]
        w_ref = lib.wav
    if wavlim is not None:
        w_ref, s_ref, serr_ref = specmatchio.truncate_spectrum(wavlim, w_ref, s_ref, serr_ref)
    
    # find unshifted spectrum
    param, spectrum = lib[lib_index]
    unshifted_file = UNSHIFTED_PATH.format(param.lib_obs)
    
    if offset:
        plt.plot(w_ref, s_ref, label="Reference")
        plot_hires_spectrum(unshifted_file, wavlim, label="Unshifted", offset=1)
        plot_library_spectrum(lib, param.lib_obs, wavlim, label="Shifted", offset=-1)
    else:
        plt.plot(w_ref, s_ref, label="Reference")
        plot_hires_spectrum(unshifted_file, wavlim, label="Unshifted")
        plot_library_spectrum(lib, param.lib_obs, wavlim, label="Shifted")
    
    plt.title('Star: {0}, Spectrum: {1}\nReference: {2}'.format(param.cps_name, param.lib_obs, ref))
    plt.xlabel('Wavelength (Angstroms)')
    plt.legend(loc='best')

def shift_data_plot(lib, lib_index):    
    """Plots lags for each order of the spectrum
    """
    param, spectrum = lib[lib_index]
    shift_datafile = '/Users/samuel/SpecMatch-Emp/lib/shift_data/{0}.txt'.format(param.lib_obs)
    list_lags = np.loadtxt(shift_datafile)
    shape = np.shape(list_lags)
    list_lags = list_lags.reshape(shape[0], 3, int(shape[1]/3))
    for i in range(len(list_lags)):
        plt.plot(list_lags[i,0], list_lags[i,1], '.')
        plt.plot(list_lags[i,0], list_lags[i,2], 'r-')
    # plt.title('Star: {0}, Spectrum: {1}'.format(param.cps_name, param.lib_obs))
    plt.xlabel('Pixel value')
    plt.ylabel('Lag')

def shift_plot(lib, lib_index, wavlim=(5160,5190)):
    gs = gridspec.GridSpec(3,1)
    plt.subplot(gs[0])
    ref_path = '/Users/samuel/Dropbox/SpecMatch-Emp/nso/nso_std.fits'
    shifted_spectrum_plot(lib, lib_index, wavlim, ref_path)
    plt.subplot(gs[1:2])
    shift_data_plot(lib, lib_index)

################################# Shift plots ##################################

################################# Match plots ##################################
def plot_match(mt, plot_targ=True, plot_resid=True, offset=False):
    """Plot a match object

    Args:
        mt (match.Match): Match object
        plot_targ (bool): If false, does not plot the target 
        plot_resid (bool): If true, plots the residuals between the spectrum
            and library values.
        offset (bool): If true, offsets the target, reference and modified spectra
    """
    if plot_targ:
        plt.plot(mt.w, mt.s_targ, label="Target")
    if mt.s_mod is None:
        off = 1 if offset else 0
        plt.plot(mt.w, mt.s_ref+off, label="Library")
    else:
        off = 1 if offset else 0
        plt.plot(mt.w, mt.s_mod+off, label="Modified library")
        off = 2 if offset else 0
        plt.plot(mt.w, mt.s_ref+off, label="Library")

    if plot_resid:
        plt.plot(mt.w, mt.best_residuals(), label="Residuals")

def plot_library_match(lib, targ_idx, ref_idx, plot_targ=True, plot_resid=True, offset=False):
    """Generate and plot the match object at the given indices

    Args:
        lib
        targ_idx
        ref_idx
    """
    targ_spec = lib.library_spectra[targ_idx]
    ref_spec = lib.library_spectra[ref_idx]
    mt = match.Match(lib.wav, targ_spec, ref_spec)
    mt.best_fit()
    plot_match(mt, plot_targ, plot_resid, offset)

################################# Match plots ##################################

############################# Library test plots ###############################
def plot_param_chi_squared(targ_idx, values, lib, param):
    """Plots chi-squared as a function of a given parameter

    Args:
        targ_idx (int): The library index of the target star
        values (pd.DataFrame): Dataframe containing param column and chi_squared column,
            sorted by chi_squared
        lib (library.Library): library object
        param (str): Parameter to be plotted
    """
    # sort matches by chi_squared
    values = values.sort_values(by='chi_squared')

    plt.plot(values[param], values.chi_squared,'.')
    plt.plot(values.head(10)[param], values.head(10).chi_squared, 'r.')
    plt.axvline(x=lib.library_params.loc[targ_idx][param], color='k')
    
def chi_squared_plot(targ_idx, df_match, lib, exclude_snr=0):
    """Creates a composite plot of chi-squared as a function of Teff, logg, [Fe/H]

    Args:
        targ_idx (int): The library index of the target star
        df_match (pd.DataFrame): Dataframe containing match results
        lib (library.Library): library object
        exclude_snr (float): (optional) Signal-to-noise ratio cutoff
    """
    star = lib.library_params.loc[targ_idx].cps_name
    grouped_match = df_match.groupby('targ_idx')
    cut = df_match.ix[grouped_match.groups[targ_idx]]
    cut.rename(columns={'ref_idx': 'lib_index'}, inplace=True)
    values = pd.merge(cut, lib.library_params, how='left', on='lib_index')
    
    # remove matches with poor snr
    snr_query = "snr > {0}".format(exclude_snr)
    values = values.query(snr_query)
    
    plt.suptitle('Star: {0}'.format(star))
    plt.subplot(131)
    plt.semilogy()
    plot_param_chi_squared(targ_idx, values, lib, 'Teff')
    plt.ylabel(r'$\chi^2$')
    plt.xlabel(r'$T_{eff}$ (K)')
    plt.subplot(132)
    plt.semilogy()
    plot_param_chi_squared(targ_idx, values, lib, 'logg')
    plt.xlabel(r'$\log\ g$ (dex)')
    plt.subplot(133)
    plt.semilogy()
    plot_param_chi_squared(targ_idx, values, lib, 'feh')
    plt.xlabel(r'$[Fe/H]$ (dex)')
    plt.tight_layout()

def library_comparison_plot(lib, param_x, param_y, xlabel=None, ylabel=None, ptlabels=False, suffix='_sm'):
    """Plots comparison between library and matched values.

    Args:
        lib (library.Library): library object containing SpecMatch results as param_sm
        param_x (str): Parameter to plot on x-axis
        param_y (str): Parameter to plot on y-axis
        xlabel (str): (optional) x-axis label
        ylabel (str): (optional) y-axis label
        ptlabels (bool): (optional) Set to true to print the name of the star next to each point
    """ 
    plt.plot(lib.library_params[param_x], lib.library_params[param_y], 'ko', label='Library value')
    x = lib.library_params[[param_x+suffix, param_x]]
    y = lib.library_params[[param_y+suffix, param_y]]
    plt.plot(x.T, y.T, 'r')
    plt.plot(x.iloc[0], y.iloc[0], 'r', label='SpecMatch-Emp value')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(loc='best')

    if ptlabels:
        lib.library_params.apply(lambda x : plt.text(x[param_x],x[param_y],x['lib_index'], size='x-small', zorder=0),  axis=1)

def library_difference_plot(lib, param, label=None, clipping=None, suffix='_sm'):
    resid = lib.library_params[param+suffix] - lib.library_params[param]
    sig = np.std(resid)
    if clipping is None:
        mask = np.full_like(resid, True, dtype=bool)
    else:
        mask = (resid < clipping*sig) & (resid > -clipping*sig)
    
    plt.plot(lib.library_params[param][mask], resid[mask], 'bo')

    mean = np.mean(resid[mask])
    rms = np.sqrt(np.mean(resid[mask]**2))
    
    ax = plt.gca()
    if clipping is None:
        plt.text(0.05, 0.1, "Mean Diff: {0:.3g}\nRMS Diff: {1:.3g}".format(mean, rms)\
            ,transform=ax.transAxes)
    else:
        plt.text(0.05, 0.1, "Mean Diff: {0:.3g}\nRMS Diff: {1:.3g}\nClipping: {2:d}".format(mean, rms, clipping)\
            +r'$\sigma$',transform=ax.transAxes)
    plt.axhline(y=0, color='k', linestyle='dashed')

    if label is not None:
        plt.xlabel(label)
        plt.ylabel(r'$\Delta\ $'+label)

def diagnostic_plots(lib, query=None, clipping=2, suffix='_sm'):
    temp_params = lib.library_params
    if query is not None:
        lib.library_params = lib.library_params.query(query)

    gs = gridspec.GridSpec(6,2)
    ax = plt.subplot(gs[0:3,0])
    library_comparison_plot(lib, 'Teff', 'logr', r'$T_{eff}$ (K)', r'$\log\ R (R_\odot)$', suffix=suffix)
    plt.xlim(plt.xlim()[::-1])
    # ax.set_yscale('log')
    ax = plt.subplot(gs[3:6,0])
    library_comparison_plot(lib, 'feh', 'logr', r'$[Fe/H]$ (dex)', r'$\log\ R (R_\odot)$', suffix=suffix)
    # ax.set_yscale('log')
    plt.subplot(gs[0:2,1])
    library_difference_plot(lib, 'Teff', r'$T_{eff}$ (K)', clipping=clipping, suffix=suffix)
    plt.xlim(plt.xlim()[::-1])
    ax = plt.subplot(gs[2:4,1])
    library_difference_plot(lib, 'logr', r'$\log\ R (R_\odot)$', clipping=clipping, suffix=suffix)
    # ax.set_yscale('log')
    plt.subplot(gs[4:6,1])
    library_difference_plot(lib, 'feh', r'$[Fe/H]$ (dex)', clipping=clipping, suffix=suffix)

    lib.library_params = temp_params


############################# Library test plots ###############################
