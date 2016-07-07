"""
@filename plotting.py

Helper functions to plot various data from SpecMatch-Emp
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import lmfit

from specmatchemp import library
from specmatchemp import match
from specmatchemp.io import specmatchio

######################### Library and Spectrum plots ###########################
def plot_library_params(lib, param_x, param_y, grouped=True):
    """Plot a H-R diagram from the library

    Args:
        lib (library.Library): The library object
        param_x (str): Parameter to be plot on the x-axis
        param_y (str): Parameter to be plot on the y-axis
        grouped (bool): (optional): Whether to group library by source catalog
    """
    x = param_x
    y = param_y
    assert x in lib.library_params.columns, "{0} not in library_params".format(x)
    assert y in lib.library_params.columns, "{0} not in library_params".format(y)
    
    if grouped:
        g = lib.library_params.groupby('source')
        for source in g.groups:
            cut = lib.library_params.ix[g.groups[source]]
            plt.plot(cut[x], cut[y], '.', label=source)
    else:
        plt.plot(lib.library_params[x], lib.library_params[y], '.')

def plot_hires_spectrum(filename, wavlim, label=None, offset=0):
    """Plot a HIRES spectrum within a given wavelength range

    Args:
        filename (str): FITS file containing HIRES spectrum
        wavlim (2-element iterable): Wavelength range
        label (str): (optional) Label for this specturm
        offset (float): (optional) Vertical offset of spectrum
    """
    w, s, serr, hdr = specmatchio.read_hires_spectrum(filename)
    w = w.reshape(-1)
    s = s.reshape(-1)
    w, s, serr = specmatchio.truncate_spectrum(wavlim, w, s)
    percen = np.percentile(s,95)
    s /= percen
    plt.plot(w, s+offset, label=label)
    plt.xlim(wavlim)

def plot_standard_spectrum(filename, wavlim, label=None, offset=0):
    """Plot a standard spectrum saved by specmatchemp.io.io within
    a given wavelength range

    Args:
        filename (str): FITS file containing spectrum
        wavlim (2-element iterable): Wavelength range
        label (str): (optional) Label for this specturm
        offset (float): (optional) Vertical offset of spectrum
    """
    w, s, serr, hdr = specmatchio.read_standard_spectrum(filename)
    w, s, serr = specmatchio.truncate_spectrum(wavlim, w, s)
    percen = np.percentile(s,95)
    s /= percen
    plt.plot(w, s+offset, label=label)
    plt.xlim(wavlim)

def plot_library_spectrum(lib, lib_index, wavlim, label=None, offset=0):
    """Plot a spectrum from the library.

    Args:
        lib (library.Library): The library object
        lib_index (int): Index of spectrum in library
        wavlim (2-element iterable): Wavelength range
        label (str): (optional) Label for this specturm
        offset (float): (optional) Vertical offset of spectrum
    """
    if lib_index not in lib:
        print("Index {0} not in library.".format(lib_index))
        return

    s = lib.library_spectra[lib_index, 0]
    serr = lib.library_spectra[lib_index,1]
    w = lib.wav
    w, s, serr = specmatchio.truncate_spectrum(wavlim, w, s, serr)
    plt.plot(w, s+offset, label=label)
    plt.xlim(wavlim)
######################### Library and Spectrum plots ###########################


################################# Shift plots ##################################
def plot_shifted_spectrum(lib, lib_index, ref, wavlim, offset=True):
    """Plot the shifted and unshifted spectra against the reference

    Args:
        lib (library.Library): The library object
        lib_index (int): Index of spectrum in library
        ref: Either path to a FITS file containing a standard spectrum, or
            a library index to the spectrum used as a reference.
        wavlim (2-element iterable): Wavelength range
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
        w_ref, s_ref, serr_ref = specmatchio.truncate_spectrum(wavlim, w_ref, s_ref, serr_ref)
    
    # find unshifted spectrum
    param, spectrum = lib[lib_index]
    unshifted_file = '/Users/samuel/Dropbox/SpecMatch-Emp/spectra/iodfitsdb/{0}.fits'.format(param.lib_obs)
    
    if offset:
        plt.plot(w_ref, s_ref, label="Reference")
        plot_hires_spectrum(unshifted_file, wavlim, label="Unshifted", offset=1)
        plot_library_spectrum(lib, param.lib_obs, wavlim, label="Shifted", offset=-1)
    else:
        plt.plot(w_ref, s_ref, label="Reference")
        plot_hires_spectrum(unshifted_file, wavlim, label="Unshifted")
        plot_library_spectrum(lib, param.lib_obs, wavlim, label="Shifted")
    
    plt.title('Star: {0}, Spectrum: {1}, Reference: {2}'.format(param.cps_name, param.lib_obs, ref))
    plt.xlabel('Wavelength (Angstroms)')
    plt.legend()

def plot_shift_data(lib, lib_index):    
    param, spectrum = lib[lib_index]
    shift_datafile = '/Users/samuel/SpecMatch-Emp/lib/shift_data/{0}.txt'.format(param.lib_obs)
    list_lags = np.loadtxt(shift_datafile)
    shape = np.shape(list_lags)
    list_lags = list_lags.reshape(shape[0], 3, int(shape[1]/3))
    for i in range(len(list_lags)):
        plt.plot(list_lags[i,0], list_lags[i,1], '.')
        plt.plot(list_lags[i,0], list_lags[i,2], 'r-')
    plt.title('Star: {0}, Spectrum: {1}'.format(param.cps_name, param.lib_obs))
    plt.xlabel('Pixel value')
    plt.ylabel('Lag')
################################# Shift plots ##################################

################################# Match plots ##################################
def plot_match(mt, plot_resid=True, offset=False):
    """Plot a match object

    Args:
        mt (match.Match): Match object
        plot_resid (bool): If true, plots the residuals between the spectrum
            and library values.
        offset (bool): If true, offsets the target, reference and modified spectra
    """
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
    
def plot_chi_squared(targ_idx, df_match, lib, exclude_snr=0):
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
    
    plt.title('Star: {0}'.format(star))
    plt.subplot(131)
    plt.semilogy()
    plot_param_chi_squared(targ_idx, values, lib, 'Teff')
    plt.xlabel('Teff (K)')
    plt.subplot(132)
    plt.semilogy()
    plot_param_chi_squared(targ_idx, values, lib, 'logg')
    plt.xlabel('Log g')
    plt.subplot(133)
    plt.semilogy()
    plot_param_chi_squared(targ_idx, values, lib, 'feh')
    plt.xlabel('[Fe/H]')

def plot_library_comparison(df_match, lib, param_x, param_y, num_mean=1, exclude_snr=0):
    """Plots comparison between library and matched values.

    Args:
        df_match (pd.DataFrame): Dataframe containing match results
        lib (library.Library): library object
        param_x (str): Parameter to plot on x-axis
        param_y (str): Parameter to plot on y-axis
        num_mean (int): (optional) Number of best matches to average
        exclude_snr (float): (optional) Signal-to-noise-ratio cutoff
    """
    grouped_match = df_match.groupby('targ_idx')
    x = []
    y = []
    
    for targ_idx in grouped_match.groups:
        # if lib.library_params.loc[targ_idx].Teff > 4500:
        #     continue
        cut = df_match.ix[grouped_match.groups[targ_idx]]
        cut.rename(columns={'ref_idx': 'lib_index'}, inplace=True)
        values = pd.merge(cut, lib.library_params, how='left', on='lib_index')
        values = values.sort_values(by='chi_squared')
        # remove matches with poor snr
        snr_query = "snr > {0}".format(exclude_snr)
        values = values.query(snr_query)

        lib_param_x = lib.library_params.loc[targ_idx][param_x]
        matched_param_x = values.head(num_mean).mean()[param_x]
        x.append(lib_param_x)
        x.append(matched_param_x)
        x.append(None)

        lib_param_y = lib.library_params.loc[targ_idx][param_y]
        matched_param_y = values.head(num_mean).mean()[param_y]
        y.append(lib_param_y)
        y.append(matched_param_y)
        y.append(None)
        
    plt.plot(lib.library_params[param_x], lib.library_params[param_y], 'ko', label='Library value')
    # lib.library_params.apply(lambda x : plt.text(x[param_x],x[param_y],x['cps_name'], size='x-small', zorder=0),  axis=1)
    plt.plot(x, y, 'r', label='SpecMatch-Emp value')

def plot_library_differences(df_match, lib, param, num_mean=1, exclude_snr=0):
    """Plots the differences between the library and matched values for a 
    single param.

    Args:
        df_match (pd.DataFrame): Dataframe containing match results
        lib (library.Library): library object
        param_x (str): Parameter to plot on x-axis
        num_mean (int): (optional) Number of best matches to average
        exclude_snr (float): (optional) Signal-to-noise-ratio cutoff
    """
    grouped_match = df_match.groupby('targ_idx')
    x=[]
    y=[]
    
    for targ_idx in grouped_match.groups:
        # if lib.library_params.loc[targ_idx].Teff > 4500:
        #     continue
        cut = df_match.ix[grouped_match.groups[targ_idx]]
        cut.rename(columns={'ref_idx': 'lib_index'}, inplace=True)
        values = pd.merge(cut, lib.library_params, how='left', on='lib_index')
        values = values.sort_values(by='chi_squared')
        # remove matches with poor snr
        snr_query = "snr > {0}".format(exclude_snr)
        values = values.query(snr_query)

        lib_param = lib.library_params.loc[targ_idx][param]
        matched_param = values.head(num_mean).mean()[param]
        x.append(lib_param)
        y.append(matched_param-lib_param)

    x = np.array(x)
    y = np.array(y)
    # clip outliers
    # sig = np.std(y)
    # mask = np.where((y < 2*sig) & (y > -2*sig))
    # x = x[mask]
    # y = y[mask]
    rms = np.sqrt(np.mean(y**2))
    ax = plt.gca()
    plt.text(0.05, 0.1, "Mean Diff: {0:.3g}\nRMS Diff: {1:.3g}".\
        format(np.mean(y),rms), transform=ax.transAxes)
    plt.axhline(y=0,color='k',linestyle='dashed')
    plt.plot(x, y, 'bo')


def print_best_fit(targ_idx, df_match, lib, num_best=1, printed_params=None):
    """Prints the best_fit parameters for the given star

    Args:
        targ_idx (int): Library index of target
        df_match (pd.DataFrame): Dataframe containing match results
        lib (library.Library): library object
        num_best (int): (optional) Number of best matches to print
        printed_params (list): (optional) List of fit parameters to print
    """
    grouped_matched = df_match.groupby('targ_idx')
    cut = df_match.ix[grouped_matched.groups[targ_idx]]
    cut = cut.sort_values(by='chi_squared').head(num_best)
    
    param, spectrum = lib[targ_idx]
    print("Best matches for star: {0}, spectrum: {1}".format(param.cps_name, param.lib_obs))
    print("Library parameters: Teff = {0:.0f}, logg = {1:.3f}, [Fe/H] = {2:.3f}"\
         .format(param.Teff, param.logg, param.feh))
    
    i = 1
    for row in cut.itertuples():
        match_idx = row.ref_idx
        param_match, spectrum_match = lib[match_idx]
        fit_params = lmfit.Parameters()
        fit_params.loads(row.fit_params)
        print("Match {0:d}, Star: {1}, Spectrum: {2}".format(i, param_match.cps_name, param_match.lib_obs))
        print("\tMatch parameters: Teff = {0:.0f}, logg = {1:.3f}, [Fe/H] = {2:.3f}"\
             .format(param_match.Teff, param_match.logg, param_match.feh))
        if printed_params is None:
            fit_params.pretty_print()
        else:
            for p in printed_params:
                if p in fit_params:
                    print("\t{0}: {1}".format(p, fit_params[p].value))
        
        i+=1

def generate_best_match(targ_idx, df_match, lib, match_idx=1):
    """Generates match objects for the best matching spectra

    Args:
        targ_idx (int): Library index of target
        df_match (pd.DataFrame): Dataframe containing match results
        lib (library.Library): library object
        match_idx (int): The rank of the best match to print
    """
    grouped_matched = df_match.groupby('targ_idx')
    cut = df_match.ix[grouped_matched.groups[targ_idx]]
    cut = cut.sort_values(by='chi_squared')

    param_targ, spec_targ = lib[targ_idx]
    
    ref_idx = cut.iloc[match_idx].ref_idx
    param_ref, spec_ref = lib[ref_idx]

    fit_params = lmfit.Parameters()
    fit_params.loads(cut.iloc[match_idx].fit_params)

    mt = match.Match(lib.wav, spec_targ[0], spec_targ[1], spec_ref[0], spec_ref[1])
    mt.create_model(fit_params)

    return mt

############################# Library test plots ###############################
