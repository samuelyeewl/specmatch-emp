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

import re
import os
import json
from argparse import ArgumentParser

from specmatchemp import library
from specmatchemp import analysis
from specmatchemp.plotting import plots

AVG_PROPS = ['Teff', 'feh', 'radius'] 
WL_AVG = [5000, 5800]

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

    outpath = os.path.join(outdir, "lincomb{0:d}_results.pdf".format(num))
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
            plots.diagnostic_plots(lib, clipping=None, suffix=suf, trend=polycoeffs)
            fig.suptitle("SpecMatch-Emp Results " + \
                "Method: Linear Combination {0:d}".format(num) + \
                "\nWavelength: {0:d} - {1:d} A".format(wl, wl+100), fontsize=16)
            plt.tight_layout(rect=[0,0.03,1,0.95])
            pdf.savefig()
            plt.close()

            # Plot cool stars
            fig = plt.figure(figsize=(15,12))
            plots.diagnostic_plots(lib, query='Teff < 4500', clipping=None, suffix=suf, trend=polycoeffs)
            fig.suptitle(r"SpecMatch-Emp Results, $T_{eff} < 4500$ K" + \
                "Method: Linear Combination {0:d}".format(num) + \
                "\nWavelength: {0:d} - {1:d} A".format(wl, wl+100), fontsize=16)
            plt.tight_layout(rect=[0,0.03,1,0.95])
            pdf.savefig()
            plt.close()

            # Plot hot stars
            fig = plt.figure(figsize=(15,12))
            plots.diagnostic_plots(lib, query='7000 >= Teff >= 4500', clipping=None, suffix=suf, trend=polycoeffs)
            fig.suptitle(r"SpecMatch-Emp Results, $7000 \geq T_{eff} \geq 4500$ K" + \
                "Method: Linear Combination {0:d}".format(num) + \
                "\nWavelength: {0:d} - {1:d} A".format(wl, wl+100), fontsize=16)
            plt.tight_layout(rect=[0,0.03,1,0.95])
            pdf.savefig()
            plt.close()

            ### Same with linear trend removed
            fig = plt.figure(figsize=(15,12))
            plots.diagnostic_plots(lib, clipping=None, suffix=suf+'_detrend')
            fig.suptitle("SpecMatch-Emp Results " + \
                "Method: Linear Combination {0:d}, Detrended".format(num) + \
                "\nWavelength: {0:d} - {1:d} A".format(wl, wl+100), fontsize=16)
            plt.tight_layout(rect=[0,0.03,1,0.95])
            pdf.savefig()
            plt.close()

            fig = plt.figure(figsize=(15,12))
            plots.diagnostic_plots(lib, query='Teff < 4500', clipping=None, suffix=suf+'_detrend')
            fig.suptitle(r"SpecMatch-Emp Results, $T_{eff} < 4500$ K" + \
                "Method: Linear Combination {0:d}, Detrended".format(num) + \
                "\nWavelength: {0:d} - {1:d} A".format(wl, wl+100), fontsize=16)
            plt.tight_layout(rect=[0,0.03,1,0.95])
            pdf.savefig()
            plt.close()

            # Plot hot stars
            fig = plt.figure(figsize=(15,12))
            plots.diagnostic_plots(lib, query='7000 >= Teff >= 4500', clipping=None, suffix=suf+'_detrend')
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
        plots.diagnostic_plots(lib, clipping=None, suffix=suf_avg)
        fig.suptitle("SpecMatch-Emp Results " + \
            "Method: Linear Combination {0:d}".format(num) + \
            "\nAveraged over wavlengths {0:d} to {1:d} A".format(WL_AVG[0], WL_AVG[1]), fontsize=16)
        plt.tight_layout(rect=[0,0.03,1,0.95])
        pdf.savefig()
        plt.close()

        # Plot cool stars
        fig = plt.figure(figsize=(15,12))
        plots.diagnostic_plots(lib, query='Teff < 4500', clipping=None, suffix=suf_avg)
        fig.suptitle(r"SpecMatch-Emp Results, $T_{eff} < 4500$ K" + \
            "Method: Linear Combination {0:d}".format(num) + \
            "\nAveraged over wavlengths {0:d} to {1:d} A".format(WL_AVG[0], WL_AVG[1]), fontsize=16)
        plt.tight_layout(rect=[0,0.03,1,0.95])
        pdf.savefig()
        plt.close()

        # Plot hot stars
        fig = plt.figure(figsize=(15,12))
        plots.diagnostic_plots(lib, query='7000 >= Teff >= 4500', clipping=None, suffix=suf_avg)
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