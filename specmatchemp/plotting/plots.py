"""
@filename plots.py

Helper functions to plot various data from SpecMatch-Emp
"""

import numpy as np
import matplotlib
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
def plot_library_params(lib, param_x, param_y, grouped=False, ptlabels=False, plt_kw={}):
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

    if ptlabels is not False:
        params.apply(lambda x : plt.text(x[param_x],x[param_y],x[ptlabels], size='x-small', zorder=0),  axis=1)

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
def plot_shifts(s, w, s_un, w_un, s_ref, w_ref, s_nso=None, w_nso=None, wavlim=None, legend=False, labels={}):
    """Plot the shifted and unshifted spectra against the reference.

    Args:
        s, w: Shifted spectrum
        s_un, w_un: Unshifted spectrum
        s_ref, w_ref: Reference (bootstrapped) spectrum
        s_nso, w_nso: (optional) nso spectrum
        wavlim: (optional) Wavelength limits
        legend: (optional) Set to true if 
    """
    if wavlim is not None:
        w, s = specmatchio.truncate_spectrum(wavlim, w, s)
        w_un, s_un = specmatchio.truncate_spectrum(wavlim, w_un, s_un)
        w_ref, s_ref = specmatchio.truncate_spectrum(wavlim, w_ref, s_ref)
        if s_nso is not None and w_nso is not None:
            w_nso, s_nso = specmatchio.truncate_spectrum(wavlim, w_nso, s_nso)

    plt.plot(w, s, '-', label="Target (shifted)")
    plt.plot(w_un, s_un-0.6, '-', label="Target (unshifted)")
    plt.plot(w_ref, s_ref+0.3, '-', label="Reference")
    if s_nso is not None and w_nso is not None:
        plt.plot(w_nso, s_nso+0.6, '-', label="NSO")

    xlim = plt.xlim()
    if 'targ_label' in labels:
        plt.text(xlim[0]+0.1, 1.05, 'Target (shifted): {0}'.format(labels['targ_label']))
        plt.text(xlim[0]+0.1, 0.45, 'Target (unshifted): {0}'.format(labels['targ_label']))
    if 'ref_label' in labels:
        plt.text(xlim[0]+0.1, 1.35, 'Reference: {0}'.format(labels['ref_label']))
    if 'nso_label' in labels:
        plt.text(xlim[0]+0.1, 1.65, 'NSO')

    plt.ylim(-0.5,1.9)
    ax = plt.gca()
    ax.axes.get_yaxis().set_ticks([])


    if legend:
        plt.legend(loc='lower left')
    plt.xlabel('Wavelength (Angstroms)')

def plot_lags(lags, center_pix, fits, legend=True):
    num_orders = lags.shape[0]
    # set different colors for each set
    colormap = plt.cm.nipy_spectral

    for i in range(num_orders):
        plt.plot(center_pix[i], lags[i], 'o', color=colormap(0.9*i/num_orders))
        plt.plot(center_pix[i], fits[i], '-', color=colormap(0.9*i/num_orders), label='{0:d}'.format(i+1))

    plt.xlabel('Pixel number')
    plt.ylabel('Shift (pixels)')
    plt.legend(loc='upper left', ncol=2, fontsize='small')


################################# Shift plots ##################################

################################# Match plots ##################################
def plot_match(mt, plot_targ=True, plot_resid=True, offset=True, labels={}, plt_kw={}):
    """Plot a match object

    Args:
        mt (match.Match): Match object
        plot_targ (bool): If false, does not plot the target 
        plot_resid (bool): If true, plots the residuals between the spectrum
            and library values.
        offset (bool): If true, offsets the target, reference and modified spectra
    """
    ax = plt.gca()
    xlim = ax.get_xlim()
    if plot_targ:
        plt.plot(mt.w, mt.s_targ, label="Target", **plt_kw)
        if 'targ_label' in labels:
            plt.text(xlim[0]+0.1, 1.05, 'Target: '+labels['targ_label'], fontsize='small')
        else:
            plt.text(xlim[0]+0.1, 1.05, 'Target', fontsize='small')
    
    if mt.s_mod is None:
        off = 1 if offset else 0
        plt.plot(mt.w, mt.s_ref+off, label="Library", **plt_kw)
        if 'ref_label' in labels:
            plt.text(xlim[0]+0.1, 1.05, 'Reference: '.format(labels['ref_label']), fontsize='small')
        else:
            plt.text(xlim[0]+0.1, 1.05, 'Reference', fontsize='small')

    else:
        off = 1 if offset else 0
        plt.plot(mt.w, mt.s_mod+off, label="Modified library", **plt_kw)
        if 'mod_label' in labels:
            plt.text(xlim[0]+0.1, 2.05, 'Reference (Modified): '+labels['mod_label'], fontsize='small')
        else:
            plt.text(xlim[0]+0.1, 2.05, 'Reference (Modified)', fontsize='small')

        off = 2 if offset else 0
        plt.plot(mt.w, mt.s_ref+off, label="Library", **plt_kw)
        if 'ref_label' in labels:
            plt.text(xlim[0]+0.1, 3.05, 'Reference: '+labels['ref_label'], fontsize='small')
        else:
            plt.text(xlim[0]+0.1, 3.05, 'Reference', fontsize='small')

    if plot_resid:
        plt.plot(mt.w, mt.best_residuals(), label="Residuals")
        if 'res_label' in labels:
            plt.text(xlim[0]+0.1, 0.05, 'Residuals: '+labels['res_label'], fontsize='small')
        else:
            plt.text(xlim[0]+0.1, 0.05, 'Residuals', fontsize='small')

    # Hide y label ticks
    ax = plt.gca()
    ax.axes.get_yaxis().set_ticks([])

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

def bestmatch_comparison_plot(res, paramx, paramy, num_best, cscol, distcol='dist'):
    res = res.sort_values(by=distcol)
    plt.plot(res[paramx], res[paramy], '.', label='_nolegend_')
    plt.plot(res.iloc[0][paramx], res.iloc[0][paramy], '*', label='Target', ms=15)
    plt.plot(res.iloc[1:num_best+1][paramx], res.iloc[1:num_best+1][paramy], 's', label='Closest stars')
    res = res.sort_values(by=cscol)
    plt.plot(res.iloc[0:num_best][paramx], res.iloc[0:num_best][paramy], '^', label='Best matches')
    plt.legend(numpoints=1, fontsize='small', loc='upper left')


############################# Library test plots ###############################
# def plot_param_chi_squared(targ_idx, values, lib, param):
#     """Plots chi-squared as a function of a given parameter

#     Args:
#         targ_idx (int): The library index of the target star
#         values (pd.DataFrame): Dataframe containing param column and chi_squared column,
#             sorted by chi_squared
#         lib (library.Library): library object
#         param (str): Parameter to be plotted
#     """
#     # sort matches by chi_squared
#     values = values.sort_values(by='chi_squared_5100')

#     plt.plot(values[param], values.chi_squared,'.')
#     plt.plot(values.head(10)[param], values.head(10).chi_squared, 'r.')
#     plt.axvline(x=lib.library_params.loc[targ_idx][param], color='k')

def plot_param_chi_squared(res, param, targ_idx, suffix):
    cs_col = 'chi_squared'+suffix
    res = res.sort_values(by=cs_col)

    plt.plot(res[param], res[cs_col], '.')
    plt.plot(res.head(10)[param], res.head(10)[cs_col], 'r.')
    plt.axvline(x=res.loc[targ_idx, param], color='k')
    

def chi_squared_plot(res, targ_idx, suffix):
    plt.subplot(131)
    plt.semilogy()
    plot_param_chi_squared(res, 'Teff', targ_idx, suffix)
    plt.ylabel(r'$\chi^2$')
    plt.xlabel(r'$T_{eff}$ (K)')
    plt.subplot(132)
    plt.semilogy()
    plot_param_chi_squared(res, 'radius', targ_idx, suffix)
    ax = plt.gca()
    ax.set_xscale('log')
    plt.xlabel(r'$R\ (R_\odot)$')
    plt.subplot(133)
    plt.semilogy()
    plot_param_chi_squared(res, 'feh', targ_idx, suffix)
    plt.xlabel(r'$[Fe/H]$ (dex)')
    

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
    resid = lib.library_params[param+suffix+'_resid']

    # sigma clipping
    sig = np.std(resid)
    if clipping is None:
        mask = np.full_like(resid, True, dtype=bool)
    else:
        mask = (resid < clipping*sig) & (resid > -clipping*sig)
    
    if param == 'dr_r':
        plt.plot(lib.library_params['radius'][mask], resid[mask], 'bo')
    else:
        plt.plot(lib.library_params[param][mask], resid[mask], 'bo')


    mean = np.mean(resid[mask])
    mean = 0 if np.isclose(mean, 0, atol=1e-4) else mean
    rms = np.sqrt(np.mean(resid[mask]**2))
    
    ax = plt.gca()
    if clipping is None:
        plt.text(0.05, 0.1, "Mean Diff: {0:.3g}\nRMS Diff: {1:.3g}".format(mean, rms)\
            ,transform=ax.transAxes)
    else:
        plt.text(0.05, 0.1, "Mean Diff: {0:.3g}\nRMS Diff: {1:.3g}\nClipping: {2:d}".format(mean, rms, clipping)\
            +r'$\sigma$',transform=ax.transAxes)
    plt.axhline(y=0, color='k', linestyle='dashed')

    if param == 'dr_r':
        plt.xlabel(label)
        plt.ylabel(r'$\Delta R/R$')
    elif label is not None:
        plt.xlabel(label)
        plt.ylabel(r'$\Delta\ $'+label)

def diagnostic_plots(lib, query=None, clipping=None, suffix='_sm', trend=None):
    temp_params = lib.library_params
    if query is not None:
        lib.library_params = lib.library_params.query(query)

    gs = gridspec.GridSpec(6,2)
    ## HR diagram
    ax = plt.subplot(gs[0:3,0])
    library_comparison_plot(lib, 'Teff', 'radius', r'$T_{eff}$ (K)', r'$R\ (R_\odot)$', suffix=suffix)
    reverse_x()
    ax.set_yscale('log')

    ax = plt.subplot(gs[3:6,0])
    library_comparison_plot(lib, 'feh', 'radius', r'$[Fe/H]$ (dex)', r'$R\ (R_\odot)$', suffix=suffix)
    ax.set_yscale('log')
    
    ## Difference plots
    plt.subplot(gs[0:2,1])
    library_difference_plot(lib, 'Teff', r'$T_{eff}$ (K)', clipping=clipping, suffix=suffix)
    # Plot trend
    if trend is not None:
        # piecewise linear trend in Teff
        # hot stars
        hot = lib.library_params.query('Teff >= 4500')
        if len(hot) >= 2:
            p = trend['Teff_hot']
            xpts = np.array([7500, 4500])
            plt.plot(xpts, p[0]*xpts+p[1], 'r-')
            ax = plt.gca()
            plt.text(0.05, 0.9, '{0:.3g}x + {1:.3g}'.format(p[0], p[1]), transform=ax.transAxes)
        
        # cool stars
        cool = lib.library_params.query('Teff <4500')
        if len(cool) >= 2:
            p = trend['Teff_cool']
            xpts = np.array([4500, 3050])
            plt.plot(xpts, p[0]*xpts+p[1], 'r-')
            ax = plt.gca()
            plt.text(0.75, 0.9, '{0:.3g}x + {1:.3g}'.format(p[0], p[1]), transform=ax.transAxes)

    reverse_x()
    
    ax = plt.subplot(gs[2:4,1])
    library_difference_plot(lib, 'dr_r', r'$R (R_\odot)$', clipping=clipping, suffix=suffix)
    if trend is not None:
        # linear trend in radius for giants
        giants = lib.library_params.query('1. < radius < 2.5')
        if len(giants) >= 2:
            p = trend['radius_giants']
            xpts = np.array([1.1, 2.5])
            plt.plot(xpts, p[0]*np.log(xpts)+p[1], 'r-')
            ax = plt.gca()
            plt.text(0.05, 0.9, '{0:.3g}x + {1:3g}'.format(p[0], p[1]), transform=ax.transAxes)

    ax.set_xscale('log')
    
    plt.subplot(gs[4:6,1])
    library_difference_plot(lib, 'feh', r'$[Fe/H]$ (dex)', clipping=clipping, suffix=suffix)
    if trend is not None:
        # linear trend in feh
        p = trend['feh']
        xpts = np.array([-0.5,0.5])
        plt.plot(xpts, p[0]*xpts+p[1], 'r-')
        ax = plt.gca()
        plt.text(0.05, 0.9, '{0:.3g}x + {1:.3g}'.format(p[0], p[1]), transform=ax.transAxes)

    lib.library_params = temp_params


############################# Library test plots ###############################
