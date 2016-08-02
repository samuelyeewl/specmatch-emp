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

import re
import os
from argparse import ArgumentParser

from specmatchemp import library
from specmatchemp import analysis
from specmatchemp.plotting import plots

def main(libpath, respath, outpath, suffix=""):
    lib = library.read_hdf(libpath, wavlim='none')
    lib.library_params.loc[:,'closest_star'] = lib.library_params.apply(analysis.find_closest_star, args=(lib,), axis=1)
    lib.library_params.loc[:,'found_closest'] = False
    res_global = pd.read_csv(respath, index_col=0)

    cscols = [c for c in list(res_global.columns) if re.search('chi_squared_\d+', c)]
    wls = [re.search('chi_squared_(\d+)$', c).group(1) for c in cscols]
    wls = map(int, wls)

    with PdfPages(outpath) as pdf:
        for (wl, cs) in zip(wls, cscols):
            # generate best match values
            suf = '_bm{0:d}'.format(wl)
            lib.library_params = analysis.generate_sm_values(lib.library_params, res_global,\
                method='best_match', suffix=suf, cscol=cs)
            lib.library_params = analysis.generate_residuals(lib.library_params, suffix=suf)
            lib.library_params.loc[:,'found_closest_{0:d}'.format(wl)] = \
                lib.library_params.closest_star == lib.library_params['best_match'+suf]
            lib.library_params.loc[:,'found_closest'] |= lib.library_params.loc[:,'found_closest_{0:d}'.format(wl)]

            # Plot entire library
            fig = plt.figure(figsize=(15,12))
            plots.diagnostic_plots(lib, clipping=None, suffix=suf)
            fig.suptitle("SpecMatch-Emp Results, Method: Best Match"+\
                "\nWavelength: {0:d} - {1:d} A".format(wl, wl+100), fontsize=16)
            plt.tight_layout(rect=[0,0.03,1,0.95])
            pdf.savefig()
            plt.close()

            # Plot cool stars
            fig = plt.figure(figsize=(15,12))
            plots.diagnostic_plots(lib, query='Teff < 4500', clipping=None, suffix=suf)
            fig.suptitle(r"SpecMatch-Emp Results, $T_{eff} < 4500$ K, Method: Best Match"+\
                "\nWavelength: {0:d} - {1:d} A".format(wl, wl+100), fontsize=16)
            plt.tight_layout(rect=[0,0.03,1,0.95])
            pdf.savefig()
            plt.close()

            # Plot hot stars
            fig = plt.figure(figsize=(15,12))
            plots.diagnostic_plots(lib, query='7000 >= Teff >= 4500', clipping=None, suffix=suf)
            fig.suptitle(r"SpecMatch-Emp Results, $7000 \geq T_{eff} \geq 4500$ K, Method: Best Match"+\
                "\nWavelength: {0:d} - {1:d} A".format(wl, wl+100), fontsize=16)
            plt.tight_layout(rect=[0,0.03,1,0.95])
            pdf.savefig()
            plt.close()

            # High stars for which closest match was found
            fig = plt.figure(figsize=(12,8))
            plots.library_comparison_plot(lib, 'Teff', 'radius', r'$T_{eff}$ (K)', r'$R\ (R_\odot)$', suffix=suf)
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
        plots.plot_library_params(lib, 'Teff', 'radius', grouped=False)
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
    psr.add_argument('outpath', type=str, help="Path to output directory")
    psr.add_argument('-s', '--suffix', type=str, default="", help="Suffix on plot files")
    args = psr.parse_args()

    main(args.library, args.respath, args.outpath, args.suffix)