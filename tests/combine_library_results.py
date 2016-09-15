#!/usr/bin/env python
"""
@filename combine_library_results.py

Combine results from internal cross-validation of library
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import os
import re
from argparse import ArgumentParser

from specmatchemp import library
from specmatchemp import specmatch
from specmatchemp import diagplots


def main(resdir, suffix, plots):
    lib = library.read_hdf(wavlim='none')

    # Create global dataframe to store results
    match_res_global = pd.DataFrame()
    lincomb_res_global = pd.DataFrame()

    for idx, row in lib.library_params.iterrows():
        targ_name = row['cps_name']
        res_path = os.path.join(resdir, targ_name + '/' + targ_name +
                                args.suffix + '_lincomb_sm.hdf')
        if not os.path.exists(res_path):
            res_path = os.path.join(resdir, targ_name + '/' + targ_name +
                                    args.suffix + '_sm.hdf')
            if not os.path.exists(res_path):
                print("Results file for {0} not found!".format(targ_name))
                continue

        print("Reading in results for {0}".format(targ_name))

        sm = specmatch.SpecMatch.read_hdf(res_path, lib)

        # Read in library grid search results
        sm.match_results['targ_idx'] = idx
        sm.match_results['ref_idx'] = sm.match_results.index
        # Concatenate library grid search results
        if match_res_global.empty:
            match_res_global = sm.match_results
        else:
            match_res_global = pd.concat((match_res_global, sm.match_results),
                                         ignore_index=True)

        if lincomb_res_global.empty:
            # Create global results dataframe containing final lincomb
            # chi-squared and generated results
            lincomb_res_global = lib.library_params.copy()
            for i in range(len(sm.lincomb_regions)):
                reg = sm.lincomb_regions[i]
                suf = '_lincomb_{0:.0f}'.format(reg[0])
                cs_col = 'chi_squared' + suf
                lincomb_res_global[cs_col] = np.nan
                for p in library.Library.STAR_PROPS:
                    p_col = p + suf
                    lincomb_res_global[p_col] = np.nan

            suf = '_lincomb'
            for p in library.Library.STAR_PROPS:
                p_col = p + suf
                p_nodetrend_col = p + '_nodetrend' + suf
                lincomb_res_global[p_col] = np.nan
                lincomb_res_global['u_'+p_col] = np.nan
                lincomb_res_global[p_nodetrend_col] = np.nan

        # Read in lincomb results for each wavelength region
        for i in range(len(sm.lincomb_regions)):
            mt = sm.lincomb_matches[i]
            reg = sm.lincomb_regions[i]
            suf = '_lincomb_{0:.0f}'.format(reg[0])
            cs_col = 'chi_squared' + suf
            lincomb_res_global.loc[idx, cs_col] = mt.best_chisq
            for p in library.Library.STAR_PROPS:
                p_col = p + suf
                lincomb_res_global.loc[idx, p_col] = sm.lincomb_results[i][p]

        # Read in averaged results
        suf = '_lincomb'
        for p in library.Library.STAR_PROPS:
            p_col = p + suf
            p_nodetrend_col = p + '_nodetrend' + suf
            lincomb_res_global.loc[idx, p_col] = sm.results[p]
            lincomb_res_global.loc[idx, 'u_'+p_col] = sm.results['u_'+p]
            lincomb_res_global.loc[idx, p_nodetrend_col] = \
                sm.results_nodetrend[p]

    # Save results
    outpath = os.path.join(resdir, 'match_results' + suffix + '.csv')
    match_res_global.index.names = ['idx']
    match_res_global.to_csv(outpath)
    outpath = os.path.join(resdir, 'lincomb_results' + suffix + '.csv')
    lincomb_res_global.to_csv(outpath)

    if plots:
        # Plots for individual regions
        plotspath = os.path.join(resdir, 'library_plots_regions.pdf')
        with PdfPages(plotspath) as pdf:
            cols = list(lincomb_res_global.columns)
            cscols = cscols = [c for c in cols if re.search(
                'chi_squared_lincomb_\d+', c)]
            regions = list(set([re.search('chi_squared_lincomb_(\d+)$', c)
                               .group(1) for c in cscols]))

            for reg in regions:
                reg = float(reg)
                suf = '_lincomb_{0:.0f}'.format(reg)
                title = 'Wavelength region: {0:.0f} - {1:.0f} A'\
                        .format(reg, reg + 100)
                plot_diag(lincomb_res_global, suf, pdf, title)

        # Plots for nodetrend parameters
        plotspath = os.path.join(resdir, 'library_plots_nodetrend.pdf')
        with PdfPages(plotspath) as pdf:
            suf = '_nodetrend_lincomb'
            title = 'Averaged over all wavelengths, not detrended'
            # Plot all stars
            plot_diag(lincomb_res_global, suf, pdf, title, trend=True)

        # Plots for detrended parameters
        plotspath = os.path.join(resdir, 'library_plots_lincomb.pdf')
        with PdfPages(plotspath) as pdf:
            suf = '_lincomb'
            title = 'Averaged over all wavelengths, detrended'
            # Plot all stars
            plot_diag(lincomb_res_global, suf, pdf, title)


def plot_diag(results, suffix, pdf, title="", trend=False):
    # Plot all stars
    fig = plt.figure(figsize=(15, 12))
    diagplots.five_pane(results, suffix, trend=trend)
    plt.suptitle('Results from Lincomb Approach\n' + title)
    plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])
    pdf.savefig()
    plt.close()

    # Plot cool stars
    fig = plt.figure(figsize=(15, 12))
    diagplots.five_pane(results.query('Teff < 4500'), suffix, trend=trend)
    plt.suptitle('Results from Lincomb Approach, Teff < 4500\n' + title)
    plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])
    pdf.savefig()
    plt.close()

    # Plot cool stars
    fig = plt.figure(figsize=(15, 12))
    diagplots.five_pane(results.query('Teff >= 4500'), suffix, trend=trend)
    plt.suptitle('Results from Lincomb Approach, Teff >= 4500\n' + title)
    plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])
    pdf.savefig()
    plt.close()


if __name__ == '__main__':
    # Argument parser
    psr = ArgumentParser(description="Combine library results")
    psr.add_argument('resdir', type=str, help="Path to results directory")
    psr.add_argument('-s', '--suffix', type=str, default="",
                     help="Suffix on individual results files")
    psr.add_argument('-p', '--plots', action='store_true',
                     help="Whether to create diagnostic plots")
    args = psr.parse_args()

    main(args.resdir, args.suffix, args.plots)
