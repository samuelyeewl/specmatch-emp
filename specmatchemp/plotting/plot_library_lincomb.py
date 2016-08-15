#!/usr/bin/env python
"""
@filename plot_library_lincomb.h5

Plots global lincomb match plots
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
import json
from argparse import ArgumentParser

from specmatchemp import library
from specmatchemp import analysis
from specmatchemp import plots

AVG_PROPS = ['Teff', 'feh', 'radius'] 
WL_AVG = [5000, 5800]

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


def main(libpath, respath, outdir, num, suffix=""):
    lib = library.read_hdf(libpath, wavlim='none')
    res_global = pd.read_csv(respath, index_col=0)

    # get wavelength regions and lincomb numbers used
    cscols = [c for c in list(res_global.columns) if re.search('chi_squared_lincomb{0:d}_\d+'.format(num), c)]
    wls = list(set([re.search('chi_squared_lincomb{0:d}_(\d+)$'.format(num), c).group(1) for c in cscols]))
    wls.sort()
    wls = list(map(int, wls))

    # number of wavelength regions being averaged over
    wls_avg = [wl for wl in wls if wl > WL_AVG[0] and wl < WL_AVG[1]]
    num_wls = len(wls_avg)

    suf_avg = '_lc{0:d}_avg'.format(num)
    for p in AVG_PROPS:
        pavg = p+suf_avg
        lib.library_params[pavg] = 0

    outpath = os.path.join(outdir, "lincomb{0:d}_results{1}.pdf".format(num, suffix))
    with PdfPages(outpath) as pdf:
        for wl in wls:
            refcol = 'ref_idxs_lincomb{0:d}_{1:d}'.format(num, wl)
            coeffcol = 'coeffs_lincomb{0:d}_{1:d}'.format(num, wl)
            # convert string to float arrays
            res_global[refcol] = res_global[refcol].apply(json.loads)
            res_global[coeffcol] = res_global[coeffcol].apply(json.loads)
            suf = '_lc{0:d}_{1:d}'.format(num, wl)

            lib.library_params = analysis.generate_sm_values(lib.library_params, res_global, \
                method='lincomb', suffix=suf, refcol=refcol, coeffcol=coeffcol)
            lib.library_params = analysis.generate_residuals(lib.library_params, suf)
            lib.library_params, polycoeffs = analysis.detrend_params(lib.library_params, suf)

            if wl in wls_avg:
                for p in AVG_PROPS:
                    pavg = p+suf_avg
                    pwl = p+suf+'_detrend'
                    lib.library_params.loc[:,pavg] += lib.library_params[pwl]/num_wls

            # Plot entire library
            fig = plt.figure(figsize=(15,12))
            diagnostic_plots(lib, clipping=None, suffix=suf, trend=polycoeffs)
            fig.suptitle("SpecMatch-Emp Results " + \
                "Method: Linear Combination {0:d}".format(num) + \
                "\nWavelength: {0:d} - {1:d} A".format(wl, wl+100), fontsize=16)
            plt.tight_layout(rect=[0,0.03,1,0.95])
            pdf.savefig()
            plt.close()

            # Plot cool stars
            fig = plt.figure(figsize=(15,12))
            diagnostic_plots(lib, query='Teff < 4500', clipping=None, suffix=suf, trend=polycoeffs)
            fig.suptitle(r"SpecMatch-Emp Results, $T_{eff} < 4500$ K" + \
                "Method: Linear Combination {0:d}".format(num) + \
                "\nWavelength: {0:d} - {1:d} A".format(wl, wl+100), fontsize=16)
            plt.tight_layout(rect=[0,0.03,1,0.95])
            pdf.savefig()
            plt.close()

            # Plot hot stars
            fig = plt.figure(figsize=(15,12))
            diagnostic_plots(lib, query='7000 >= Teff >= 4500', clipping=None, suffix=suf, trend=polycoeffs)
            fig.suptitle(r"SpecMatch-Emp Results, $7000 \geq T_{eff} \geq 4500$ K" + \
                "Method: Linear Combination {0:d}".format(num) + \
                "\nWavelength: {0:d} - {1:d} A".format(wl, wl+100), fontsize=16)
            plt.tight_layout(rect=[0,0.03,1,0.95])
            pdf.savefig()
            plt.close()

            ### Same with linear trend removed
            fig = plt.figure(figsize=(15,12))
            diagnostic_plots(lib, clipping=None, suffix=suf+'_detrend')
            fig.suptitle("SpecMatch-Emp Results " + \
                "Method: Linear Combination {0:d}, Detrended".format(num) + \
                "\nWavelength: {0:d} - {1:d} A".format(wl, wl+100), fontsize=16)
            plt.tight_layout(rect=[0,0.03,1,0.95])
            pdf.savefig()
            plt.close()

            fig = plt.figure(figsize=(15,12))
            diagnostic_plots(lib, query='Teff < 4500', clipping=None, suffix=suf+'_detrend')
            fig.suptitle(r"SpecMatch-Emp Results, $T_{eff} < 4500$ K" + \
                "Method: Linear Combination {0:d}, Detrended".format(num) + \
                "\nWavelength: {0:d} - {1:d} A".format(wl, wl+100), fontsize=16)
            plt.tight_layout(rect=[0,0.03,1,0.95])
            pdf.savefig()
            plt.close()

            # Plot hot stars
            fig = plt.figure(figsize=(15,12))
            diagnostic_plots(lib, query='7000 >= Teff >= 4500', clipping=None, suffix=suf+'_detrend')
            fig.suptitle(r"SpecMatch-Emp Results, $7000 \geq T_{eff} \geq 4500$ K" + \
                "Method: Linear Combination {0:d}, Detrended".format(num) + \
                "\nWavelength: {0:d} - {1:d} A".format(wl, wl+100), fontsize=16)
            plt.tight_layout(rect=[0,0.03,1,0.95])
            pdf.savefig()
            plt.close()
    
        # Same plots averaged over all wavelengths
        lib.library_params = analysis.generate_residuals(lib.library_params, suf_avg, props=AVG_PROPS)

        # Plot entire library
        fig = plt.figure(figsize=(15,12))
        diagnostic_plots(lib, clipping=None, suffix=suf_avg)
        fig.suptitle("SpecMatch-Emp Results " + \
            "Method: Linear Combination {0:d}".format(num) + \
            "\nAveraged over wavlengths {0:d} to {1:d} A".format(WL_AVG[0], WL_AVG[1]), fontsize=16)
        plt.tight_layout(rect=[0,0.03,1,0.95])
        pdf.savefig()
        plt.close()

        # Plot cool stars
        fig = plt.figure(figsize=(15,12))
        diagnostic_plots(lib, query='Teff < 4500', clipping=None, suffix=suf_avg)
        fig.suptitle(r"SpecMatch-Emp Results, $T_{eff} < 4500$ K" + \
            "Method: Linear Combination {0:d}".format(num) + \
            "\nAveraged over wavlengths {0:d} to {1:d} A".format(WL_AVG[0], WL_AVG[1]), fontsize=16)
        plt.tight_layout(rect=[0,0.03,1,0.95])
        pdf.savefig()
        plt.close()

        # Plot hot stars
        fig = plt.figure(figsize=(15,12))
        diagnostic_plots(lib, query='7000 >= Teff >= 4500', clipping=None, suffix=suf_avg)
        fig.suptitle(r"SpecMatch-Emp Results, $7000 \geq T_{eff} \geq 4500$ K" + \
            "Method: Linear Combination {0:d}".format(num) + \
            "\nAveraged over wavelengths {0:d} to {1:d} A".format(WL_AVG[0], WL_AVG[1]), fontsize=16)
        plt.tight_layout(rect=[0,0.03,1,0.95])
        pdf.savefig()
        plt.close()




if __name__ == '__main__':
    psr = ArgumentParser(description="Produce plots for analysis of match results")
    psr.add_argument('library', type=str, help="Path to library h5 file")
    psr.add_argument('respath', type=str, help="Path to results csv file")
    psr.add_argument('outdir', type=str, help="Path to output directory")
    psr.add_argument('num', type=int, help="Number of spectra used in lincomb")
    psr.add_argument('-s', '--suffix', type=str, default="", help="Suffix on plot files")
    args = psr.parse_args()

    main(args.library, args.respath, args.outdir, args.num, args.suffix)