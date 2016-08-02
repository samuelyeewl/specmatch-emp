#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd 

from specmatchemp import library
from specmatchemp.plotting import plots

import os
import sys
from argparse import ArgumentParser


def generate_sm_values(lib, results, method='best_match', cscol='chi_squared', num_mean=3, suffix='_sm'):
    params = lib.library_params
    grouped_results = results.groupby('targ_idx')

    # Use closest matched spectra as sm values
    if method == 'best_match':
        params['best_match'] = params.lib_index.apply(\
            lambda i: grouped_results.get_group(i).sort_values(by=cscol).iloc[0].ref_idx)
        for p in library.STAR_PROPS:
            psm = p+suffix
            params[psm] = params.best_match.apply(lambda i: params.loc[i, p])
        
        params['best_chi_squared'] = params.lib_index.apply(\
            lambda i: grouped_results.get_group(i).sort_values(by=cscol).iloc[0][cscol])

    elif method == 'average':
        params['best_n'] = params.lib_index.apply(lambda i: \
            np.array(grouped_results.get_group(i).sort_values(by=cscol).iloc[0:num_mean].ref_idx))
        for p in library.STAR_PROPS:
            psm = p+suffix
            params[psm] = params.best_n.apply(lambda i: params.loc[i, p].mean())
        
    return params

def dist(param1, param2):
    diff_teff = ((param1.Teff - param2.Teff)/100)**2
    diff_logg = ((param1.logg - param2.logg)/0.1)**2
    diff_feh = ((param1.feh - param2.feh)/0.1)**2
    return diff_teff + diff_logg + diff_feh

def find_closest_star(row, lib):
    return lib.library_params.apply(dist, args=(row,), axis=1).sort_values().index[1]


def main(library_path, results_path, outdir, prefix, title, minw):
    lib = library.read_hdf(library_path, wavlim='none')
    results = pd.DataFrame.from_csv(results_path, index_col=0)

    # # Plot library
    # fig = plt.figure(figsize=(12,8))
    # plots.plot_library_params(lib, 'Teff', 'logg')
    # plots.reverse_x()
    # plots.reverse_y()
    # plt.legend()
    # plt.savefig(os.path.join(outdir, prefix+"_library_Teff_logg.png"))
    # plt.close(fig)

    # fig = plt.figure(figsize=(12,8))
    # plots.plot_library_params(lib, 'feh', 'logg')
    # plots.reverse_y()
    # plt.legend()
    # plt.savefig(os.path.join(outdir, prefix+"_library_feh_logg.png"))
    # plt.close(fig)

    # Diagnostic plot for best match
    lib.library_params = generate_sm_values(lib, results, method='best_match', cscol='chi_squared_{0:d}'.format(minw), suffix='_bm')

    fig = plt.figure(figsize=(15,12))
    plots.diagnostic_plots(lib, clipping=None, suffix='_bm')
    fig.suptitle("SpecMatch-Emp Results, Method: Best Match\n"+title, fontsize=18)
    plt.tight_layout(rect=[0,0.03,1,0.95])
    plt.savefig(os.path.join(outdir, prefix+"_diagnostic.png"))
    plt.close(fig)

    fig = plt.figure(figsize=(15,12))
    plots.diagnostic_plots(lib, clipping=2, suffix='_bm')
    fig.suptitle("SpecMatch-Emp Results, Method: Best Match, Residuals clipped\n"+title, fontsize=18)
    plt.tight_layout(rect=[0,0.03,1,0.95])
    plt.savefig(os.path.join(outdir, prefix+"_diagnostic_clipped.png"))
    plt.close(fig)

    fig = plt.figure(figsize=(15,12))
    plots.diagnostic_plots(lib, query='Teff < 4500', clipping=None, suffix='_bm')
    fig.suptitle(r"SpecMatch-Emp Results, $T_{eff} < 4500$ K, Method: Best Match"+"\n"+title, fontsize=18)
    plt.tight_layout(rect=[0,0.03,1,0.95])
    plt.savefig(os.path.join(outdir, prefix+"_diagnostic_cool.png"))
    plt.close(fig)

    # Highlight points where closest star was found
    lib.library_params['closest_star'] = lib.library_params.apply(find_closest_star, args=(lib,), axis=1)
    fig = plt.figure(figsize=(12,8))
    plots.library_comparison_plot(lib, 'Teff', 'logg', xlabel=r'$T_{eff}$ (K)', ylabel=r'$\log\ g$ (dex)', suffix='_bm')
    mask = lib.library_params.closest_star == lib.library_params.best_match
    plt.plot(lib.library_params[mask].Teff,lib.library_params[mask].logg, 'o', color='cyan')
    plots.reverse_x()
    plots.reverse_y()
    plt.title("SpecMatch-Emp Results, Method: Best Match\n"+title)
    plt.savefig(os.path.join(outdir, prefix+"_best_match_highlights.png"))
    plt.close(fig)

    # Chi-squared histogram
    fig = plt.figure(figsize=(10,7))
    lib.library_params.best_chi_squared.hist()
    plt.title("Histogram of lowest chi-squared match\n"+title)
    plt.savefig(os.path.join(outdir, prefix+"_chi_squared_hist.png"))
    plt.close(fig)

    # Diagnostic plot for average
    lib.library_params = generate_sm_values(lib, results, method='average', suffix='_avg')

    fig = plt.figure(figsize=(15,12))
    plots.diagnostic_plots(lib, clipping=None, suffix='_avg')
    fig.suptitle("SpecMatch-Emp Results, Method: Average of 3\n"+title, fontsize=18)
    plt.tight_layout(rect=[0,0.03,1,0.95])
    plt.savefig(os.path.join(outdir, prefix+"_diagnostic_avg3.png"))
    plt.close(fig)

    fig = plt.figure(figsize=(15,12))
    plots.diagnostic_plots(lib, query='Teff < 4500', clipping=None, suffix='_avg')
    fig.suptitle(r"SpecMatch-Emp Results, $T_{eff} < 4500$ K, Method: Average of 3"+"\n"+title, fontsize=18)
    plt.tight_layout(rect=[0,0.03,1,0.95])
    plt.savefig(os.path.join(outdir, prefix+"_diagnostic_avg3_cool.png"))
    plt.close(fig)
    

if __name__ == '__main__':
    # Argument parser
    psr = ArgumentParser(description="Produce plots for analysis of match results")
    psr.add_argument('library', type=str, help="Path to library h5 file")
    psr.add_argument('results', type=str, help="Path to results csv file")
    psr.add_argument('outdir', type=str, help="Path to output directory")
    psr.add_argument('prefix', type=str, help="String to prefix to output files")
    psr.add_argument('title', type=str, help="String to add to title")
    psr.add_argument('minw', type=int, help="Minimum wavelength")
    args = psr.parse_args()

    mpl.rcParams['figure.dpi'] = 200


    if not os.path.isfile(args.library):
        print("Could not find {0}".format(args.library))
        sys.exit()
    if not os.path.isfile(args.results):
        print("Could not find {0}".format(args.results))
        sys.exit()

    main(args.library, args.results, args.outdir, args.prefix, args.title, args.minw)