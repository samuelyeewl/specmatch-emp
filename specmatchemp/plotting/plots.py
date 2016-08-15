"""
@filename plots.py

Helper functions to plot various data from SpecMatch-Emp
"""

from __future__ import print_function

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.transforms as transforms
import pandas as pd
import os
import lmfit


UNSHIFTED_PATH = '/Users/samuel/Dropbox/SpecMatch-Emp/spectra/iodfitsdb/{0}.fits'
def reverse_x():
    """Reverses the x-axis of the current figure"""
    plt.xlim(plt.xlim()[::-1])

def reverse_y():
    """Reverses the y-axis of the current figure"""
    plt.ylim(plt.ylim()[::-1])

def hide_x_ticks():
    """Hide x label ticks"""
    ax = plt.gca()
    ax.axes.get_xaxis().set_ticks([])

def hide_y_ticks():
    """Hide y label ticks"""
    ax = plt.gca()
    ax.axes.get_yaxis().set_ticks([])

def annotate_point(x, y, text, offset=5, offset_x=None, offset_y=None, text_kw={}):
    """Annotates the point at a given x, y position (in data coordinates),
    at a given pixel offset.

    Args:
        x: x-coordinate in data space
        y: y-coordinate in data space
        text (str): String to annotate
        offset: (optional) pixel offset to use
        offset_x, offset_y: (optional) pixel offset to use in x, y directions
        text_kw (dict): (optional) any additional keywords to pass to pyplot.text
    """
    ax = plt.gca()
    # transform to pixel coords
    disp_coords = ax.transData.transform((x, y))
    if offset_x is None or offset_y is None:
        offset_x = offset
        offset_y = offset
    disp_coords = (disp_coords[0]+offset_x, disp_coords[1]+offset_y)
    # invert transform to go back to data coords
    data_coords = ax.transData.inverted().transform(disp_coords)
    plt.text(data_coords[0], data_coords[1], text, **text_kw)

def annotate_spectrum(text, spec_offset=0, offset_x=10, offset_y=5, align='left', text_kw={}):
    """Annotates a spectrum.

    Args:
        text (str): String to annotate
        spec_offset: (optional) Vertical offset of spectrum
        offset_x: (optional) Pixel offset from left/right boundary
        offset_y: (optional) Vertical pixel offset from spectrum
        align: (optional) 'left' or 'right' alignment for text
        text_kw (dict): (optional) any additional keywords to pass to pyplot.text
    """
    ax = plt.gca()
    xlim = ax.get_xlim()
    if align == 'left':
        xpos = xlim[0]
        offset_x = abs(offset_x)
    elif align == 'right':
        xpos = xlim[1]
        offset_x = -abs(offset_x)
    else:
        return

    # transform to pixel coords
    disp_coords = ax.transData.transform((xpos, spec_offset+1))
    disp_coords = (disp_coords[0]+offset_x, disp_coords[1]+offset_y)
    # invert transform to go back to data coords
    data_coords = ax.transData.inverted().transform(disp_coords)
    ax_coords = ax.transAxes.inverted().transform(disp_coords)
    # fix y position in data coordinates (fixed offset from spectrum)
    # but allow x position to float
    trans = transforms.blended_transform_factory(ax.transAxes, ax.transData)

    bbox=dict(facecolor='white', edgecolor='none',alpha=0.8)
    plt.text(ax_coords[0], data_coords[1], text, bbox=bbox, transform=trans, horizontalalignment=align, **text_kw)
    

def set_tight_lims(data_x, data_y, center_x=None, center_y=None, mode='symmetric', buf=0.3):
    """Sets plot limits around a target subset of data, centered at
    a given point.
    
    Args:
        data_x (np.ndarray): x-coordinates of data
        data_y (np.ndarray): y-coordinates of data
        center_x (optional [float]): x-coordinate of center point
        center_y (optional [float]): y-coordinate of center point
        mode: (optional) 'symmetric': Make limits symmetric about target
                        'tight': Use asymmetric limits 
        buf (float): Buffer radius 
    """
    ax = plt.gca()

    if center_x is None:
        maxx = max(data_x)
        minx = min(data_x)
        sepx = maxx-minx
        maxx = maxx + buf*sepx
        minx = minx - buf*sepx
        ax.set_xlim((minx, maxx))
    else:
        distx = data_x - center_x
        maxx = max(max(distx), 0)
        minx = min(min(distx), 0)
        if mode == 'symmetric':
            limx = max(abs(maxx), abs(minx))
            limx = limx + buf*limx
            ax.set_xlim((center_x-limx, center_x+limx))
        elif mode == 'tight':
            maxx = maxx + buf*maxx if maxx != 0 else -buf*minx
            minx = minx + buf*minx if minx != 0 else -buf*maxx
            ax.set_xlim((center_x+minx, center_x+maxx))

    if center_y is None:
        maxy = max(data_y)
        miny = min(data_y)
        sepy = maxy-miny
        maxy = maxy + buf*sepy
        miny = miny - buf*sepy
        ax.set_ylim((miny, maxy))
    else:
        disty = data_y - center_y
        maxy = max(max(disty), 0)
        miny = min(min(disty), 0)
        if mode == 'symmetric':
            limy = max(abs(maxy), abs(miny))
            limy = limy + buf*limy
            ax.set_ylim((center_y-limy, center_y+limy))
        elif mode == 'tight':
            maxy = maxy + buf*maxy if maxy != 0 else -buf*miny
            miny = miny + buf*miny if miny != 0 else -buf*maxy
            ax.set_ylim((center_y+miny, center_y+maxy))
    

def set_lims_around_targ(library_params, paramx, paramy, targ_idx, idxs, mode='symmetric', buf=0.3):
    """Sets plot limits around a target

    Args:
        library_params (pd.DataFrame): Parameter table
        paramx, paramy (str): Parameters plot on x, y axes
        targ_idx (int): Index of target star
        idxs (iterable of ints): Index of other stars
        mode: (optional) 'symmetric': Make limits symmetric about target
                        'tight': Use asymmetric limits 
        buf (float): Buffer radius 
    """
    targx = library_params.loc[targ_idx, paramx]
    distx = library_params.loc[idxs, paramx] - targx
    ax = plt.gca()
    # get right and left limits
    maxx = max(max(distx), 0)
    minx = min(min(distx), 0)
    if mode == 'symmetric':
        limx = max(abs(maxx), abs(minx))
        limx = limx + buf*limx
        ax.set_xlim((targx-limx, targx+limx))
    elif mode == 'tight':
        maxx = maxx + buf*maxx if maxx != 0 else -buf*minx
        minx = minx + buf*minx if minx != 0 else -buf*maxx
        ax.set_xlim((targx+minx, targx+maxx))

    targy = library_params.loc[targ_idx, paramy]
    disty = library_params.loc[idxs, paramy] - targy
    ax = plt.gca()
    # get upper and lower limits
    maxy = max(max(disty), 0)
    miny = min(min(disty), 0)
    if mode == 'symmetric':
        limy = max(abs(maxy), abs(miny))
        limy = limy + buf*limy
        ax.set_ylim((targy-limy, targy+limy))
    elif mode == 'tight':
        maxy = maxy + buf*maxy if maxy != 0 else -buf*miny
        miny = miny + buf*miny if miny != 0 else -buf*maxy
        ax.set_ylim((targy+miny, targy+maxy))

def hide_y_ticks():
    # Hide y label ticks
    ax = plt.gca()
    ax.axes.get_yaxis().set_ticks([])

def label_axes(param_x, param_y):
    """Convenience function for tweaking axes to make plots
    """
    if param_x is 'Teff' and param_y is 'radius':
        yt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1, 2, 3, 4, 5, 7, 10, 20]
        plt.semilogy()
        plt.xlim(plt.xlim()[::-1])
        plt.ylim(yt[0],yt[-1])
        plt.xlabel('Effective Temperature (K)')
        plt.ylabel('Stellar Radius (Solar-radii)')
        plt.yticks(yt,yt)
        return None

    assert False, "invalid param_x, param_y"

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
    plt.plot(w_un, s_un-0.8, '-', label="Target (unshifted)")
    plt.plot(w_ref, s_ref+0.5, '-', label="Reference")
    if s_nso is not None and w_nso is not None:
        plt.plot(w_nso, s_nso+1.0, '-', label="NSO")
    resid, residw = utils.calc_residuals(s, w, s_ref, w_ref)
    plt.plot(residw, resid-1, '-', label="Residuals")

    # xlim = plt.xlim()
    # if 'targ_label' in labels:
    #     plt.text(xlim[0]+0.1, 1.05, 'Target (shifted): {0}'.format(labels['targ_label']))
    #     plt.text(xlim[0]+0.1, 0.45, 'Target (unshifted): {0}'.format(labels['targ_label']))
    # if 'ref_label' in labels:
    #     plt.text(xlim[0]+0.1, 1.35, 'Reference: {0}'.format(labels['ref_label']))
    # if 'nso_label' in labels:
    #     plt.text(xlim[0]+0.1, 1.65, 'NSO')

    plt.ylim(-1.2,2.2)
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
        plt.plot(center_pix[i], fits[i], '-', color=colormap(0.9*i/num_orders), label='{0:d}'.format(i))

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

    hide_y_ticks()

def plot_matchlincomb(mt, targ_label=None, ref_labels=None, mod_label=None, plt_kw={}):
    """Plot a MatchLincomb object

    Args:
        mt (match.MatchLincomb): MatchLincomb object
    """
    plt.plot(mt.w, mt.s_targ, label="Target", **plt_kw)
    if targ_label is not None:
        annotate_spectrum(r"Target: "+targ_label, spec_offset=0)

    for i in range(mt.num_refs):
        offset = 1.5 + i*0.5
        plt.plot(mt.w, mt.s_refs[i,0]+offset, color='0.4', label="Reference {0:d}".format(i+1))
        if ref_labels is not None:
            annotate_spectrum("Reference {0:d}: ".format(i+1)+ref_labels[i], offset)

    plt.plot(mt.w, mt.s_mod+0.5, color='r', label="Linear Combination")
    if mod_label is not None:
        annotate_spectrum("Linear Combination: "+mod_label, spec_offset=0.5)

    plt.plot(mt.w, mt.best_residuals(), color='c', label="Residuals")
    annotate_spectrum(r"Residuals", spec_offset=-1)

    plt.ylim((-0.2, 2.3+mt.num_refs*0.5))

    hide_y_ticks()


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
    plt.legend(numpoints=1, fontsize='small', loc='best')


def lincomb_refs_plot(library_params, paramx, paramy, targ_idx, ref_idxs, annot=None, zoom=False, legend=True):
    plt.plot(library_params[paramx], library_params[paramy], '.', label='_nolegend_')
    plt.plot(library_params.loc[targ_idx, paramx], library_params.loc[targ_idx, paramy], '*', label='Target', ms=15)
    plt.plot(library_params.loc[ref_idxs, paramx], library_params.loc[ref_idxs, paramy], '^', label='References')

    if zoom:
        set_lims_around_targ(library_params, paramx, paramy, targ_idx, ref_idxs)

    # add annotations
    if annot is not None:
        ax = plt.gca()
        for i, ref in enumerate(ref_idxs):
            annotate_point(library_params.loc[ref, paramx], library_params.loc[ref, paramy], annot[i], text_kw={'fontsize':'small'})

    if legend:
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

def plot_param_chi_squared(res, param, targ_idx=None, suffix=''):
    cs_col = 'chi_squared'+suffix
    res = res.sort_values(by=cs_col)

    plt.plot(res[param], res[cs_col], '.')
    plt.plot(res.head(5)[param], res.head(5)[cs_col], 'r.')
    if targ_idx is not None:
        plt.axvline(x=res.loc[targ_idx, param], color='k')


def chi_squared_plot(res, targ_idx, suffix):
    plt.subplot(131)
    plt.semilogy()
    plot_param_chi_squared(res, 'Teff', targ_idx, suffix)
    plt.ylabel(r'$\chi^2$')
    plt.xlabel(r'$T_{eff}$ (K)')
    plt.xticks([3000,4000,5000,6000,7000])
    reverse_x()
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
        lib.library_params.apply(lambda x : plt.text(x[param_x],x[param_y],' '+x['cps_name'], size='x-small', zorder=0),  axis=1)

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

def diagnostic_plots(lib, query=None, clipping=None, suffix='_sm', trend=None, ptlabels=False):
    temp_params = lib.library_params
    if query is not None:
        lib.library_params = lib.library_params.query(query)

    gs = gridspec.GridSpec(6,2)
    ## HR diagram
    ax = plt.subplot(gs[0:3,0])
    library_comparison_plot(lib, 'Teff', 'radius', r'$T_{eff}$ (K)', r'$R\ (R_\odot)$', suffix=suffix, ptlabels=ptlabels)
    reverse_x()
    ax.set_yscale('log')

    ax = plt.subplot(gs[3:6,0])
    library_comparison_plot(lib, 'feh', 'radius', r'$[Fe/H]$ (dex)', r'$R\ (R_\odot)$', suffix=suffix, ptlabels=ptlabels)
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

# from specmatchemp import library
from specmatchemp import match
from specmatchemp.io import specmatchio
from specmatchemp.utils import utils
