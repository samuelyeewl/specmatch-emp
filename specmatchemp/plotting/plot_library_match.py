#!/usr/bin/env python
"""
@filename plot_library_match.h5

Plots global match plots
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import re
import os
from argparse import ArgumentParser

from specmatchemp import library
from specmatchemp import analysis
from specmatchemp import plots

def diagnostic_plots(lib, query=None, clipping=None, suffix='_sm', trend=None, ptlabels=False):
    temp_params = lib.library_params
    if query is not None:
        lib.library_params = lib.library_params.query(query)

    gs = gridspec.GridSpec(6,2)
    ## HR diagram
    ax = plt.subplot(gs[0:3,0])
    library_comparison_plot(lib, 'Teff', 'radius', r'$T_{eff}$ (K)', r'$R\ (R_\odot)$', suffix=suffix, ptlabels=ptlabels)
    plots.reverse_x()
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

    plots.reverse_x()
    
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


def five_pane_plots(lib, res, minw, query=None, ptlabels=False):
    # generate best match values
    suf = '_bm{0:d}'.format(minw)
    cscol = 'chi_squared_{0:d}'.format(minw)
    if 'Teff'+suf not in lib.library_params.columns:
        lib.library_params = analysis.generate_sm_values(lib.library_params, res,\
                method='best_match', suffix=suf, cscol=cscol)
        lib.library_params = analysis.generate_residuals(lib.library_params, suffix=suf)

    diagnostic_plots(lib, suffix=suf, clipping=None, query=query, ptlabels=ptlabels)
    plt.suptitle("SpecMatch-Emp Results, Method: Best Match"+\
                "\nWavelength: {0:d} - {1:d} A".format(minw, minw+100), fontsize=16)
    plt.tight_layout(rect=[0,0.03,1,0.95])


def main(libpath, respath, outpath, suffix=""):
    lib = library.read_hdf(libpath, wavlim='none')
    lib.library_params.loc[:,'closest_star'] = lib.library_params.apply(analysis.find_closest_star, args=(lib,), axis=1)
    lib.library_params.loc[:,'found_closest'] = False
    res_global = pd.read_csv(respath, index_col=0)

    cscols = [c for c in list(res_global.columns) if re.search('chi_squared_\d+', c)]
    wls = [re.search('chi_squared_(\d+)$', c).group(1) for c in cscols]
    wls.sort()
    wls = map(int, wls)

    with PdfPages(outpath) as pdf:
        for (wl, cs) in zip(wls, cscols):
            # generate best match values
            suf = '_bm{0:d}'.format(wl)
            lib.library_params = analysis.generate_sm_values(lib.library_params, res_global,\
                method='best_match', suffix=suf, cscol=cs)
            lib.library_params = analysis.generate_residuals(lib.library_params, suffix=suf)
            

            # Plot entire library
            fig = plt.figure(figsize=(15,12))
            five_pane_plots(lib, res_global, wl)
            pdf.savefig()
            plt.close()

            # Plot cool stars
            fig = plt.figure(figsize=(15,12))
            five_pane_plots(lib, res_global, wl, query='Teff < 4500')
            pdf.savefig()
            plt.close()

            # Plot hot stars
            fig = plt.figure(figsize=(15,12))
            five_pane_plots(lib, res_global, wl, query='7000 >= Teff >= 4500')
            pdf.savefig()
            plt.close()

            # Highlight stars for which closest match was found
            lib.library_params.loc[:,'found_closest_{0:d}'.format(wl)] = \
                lib.library_params.closest_star == lib.library_params['best_match'+suf]
            lib.library_params.loc[:,'found_closest'] |= lib.library_params.loc[:,'found_closest_{0:d}'.format(wl)]

            fig = plt.figure(figsize=(12,8))
            library_comparison_plot(lib, 'Teff', 'radius', r'$T_{eff}$ (K)', r'$R\ (R_\odot)$', suffix=suf)
            mask = lib.library_params['found_closest_{0:d}'.format(wl)]
            plt.plot(lib.library_params[mask].Teff,lib.library_params[mask].radius, 'o', color='cyan')
            plots.reverse_x()
            ax = fig.gca()
            ax.set_yscale('log')
            plt.ylim((0.1, 16))
            plt.title("Closest star found by SpecMatch Best Match" +\
                "\nWavelength: {0:d} - {1:d} A".format(wl, wl+100), fontsize=16)
            pdf.savefig()
            plt.close()

        # Plot overall highlights
        fig = plt.figure(figsize=(12,8))
        lib.plot('Teff', 'radius')
        mask = lib.library_params['found_closest']
        plt.plot(lib.library_params[mask].Teff,lib.library_params[mask].radius, 'o', color='cyan')
        plots.reverse_x()
        ax = fig.gca()
        ax.set_yscale('log')
        plt.ylim((0.1, 16))
        plt.title("Closest star found by SpecMatch Best Match" +\
            "\nWavelength: Any", fontsize=16)
        pdf.savefig()
        plt.close()


if __name__ == '__main__':
    psr = ArgumentParser(description="Produce plots for analysis of match results")
    psr.add_argument('library', type=str, help="Path to library h5 file")
    psr.add_argument('respath', type=str, help="Path to results csv file")
    psr.add_argument('outpath', type=str, help="Path to output file")
    psr.add_argument('-s', '--suffix', type=str, default="", help="Suffix on plot files")
    args = psr.parse_args()

    main(args.library, args.respath, args.outpath, args.suffix)