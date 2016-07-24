#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd 

from specmatchemp import library
from specmatchemp import analysis
from specmatchemp.plotting import plots

import os
import sys
from argparse import ArgumentParser

def main(library_path, results_path, outdir, prefix, title, wavelengths):
    lib = library.read_hdf(library_path, 'none')
    num_wls = len(wavelengths)

    for p in library.STAR_PROPS:
        pavg = p+'_sm'
        lib.library_params[pavg] = 0

    for wl in wavelengths:
        results = pd.DataFrame.from_csv(results_path.format(wl), index_col=0)

        results['ref_idxs'] = results.ref_idxs.apply(lambda s: \
            np.array(s.strip('[]').split()).astype(int))
        results['coeffs'] = results.coeffs.apply(lambda s: \
            np.array(s.strip('[]').split()).astype(float))

        suffix = '_{0:d}_sm'.format(wl)

        lib.library_params = analysis.generate_sm_values(lib.library_params, results, method='lincomb', suffix='_{0:d}_sm'.format(wl))

        for p in library.STAR_PROPS:
            pavg = p+'_sm'
            pwl = p+'_{0:d}_sm'.format(wl)
            lib.library_params.loc[:,pavg] += lib.library_params[pwl]/num_wls

    lib.library_params = analysis.generate_residuals(lib.library_params, '_sm')
    lib.library_params, polycoeffs = analysis.detrend_params(lib.library_params, '_sm')

    # Diagnostic plots for all, cool, hot stars
    fig = plt.figure(figsize=(15,12))
    plots.diagnostic_plots(lib, clipping=None, suffix='_sm', trend=polycoeffs)
    fig.suptitle("SpecMatch-Emp Results, Method: Linear Combination\n"+title, fontsize=18)
    plt.tight_layout(rect=[0,0.03,1,0.95])
    plt.savefig(os.path.join(outdir, prefix+"_diagnostic.png"))
    plt.close(fig)

    fig = plt.figure(figsize=(15,12))
    plots.diagnostic_plots(lib, query='Teff < 4500', clipping=None, suffix='_sm', trend=polycoeffs)
    fig.suptitle(r"SpecMatch-Emp Results, $T_{eff} < 4500$ K, Method: Linear Combination"+"\n"+title, fontsize=18)
    plt.tight_layout(rect=[0,0.03,1,0.95])
    plt.savefig(os.path.join(outdir, prefix+"_diagnostic_cool.png"))
    plt.close(fig)

    fig = plt.figure(figsize=(15,12))
    plots.diagnostic_plots(lib, query='Teff > 4500', clipping=None, suffix='_sm', trend=polycoeffs)
    fig.suptitle(r"SpecMatch-Emp Results, $T_{eff} > 4500$ K, Method: Linear Combination"+"\n"+title, fontsize=18)
    plt.tight_layout(rect=[0,0.03,1,0.95])
    plt.savefig(os.path.join(outdir, prefix+"_diagnostic_hot.png"))
    plt.close(fig)

    # Diagnostic plots with linear trend removed
    fig = plt.figure(figsize=(15,12))
    plots.diagnostic_plots(lib, clipping=None, suffix='_sm_detrend')
    fig.suptitle("SpecMatch-Emp Results, Method: Linear Combination, Detrended\n"+title, fontsize=18)
    plt.tight_layout(rect=[0,0.03,1,0.95])
    plt.savefig(os.path.join(outdir, prefix+"_diagnostic_detrend.png"))
    plt.close(fig)

    fig = plt.figure(figsize=(15,12))
    plots.diagnostic_plots(lib, query='Teff < 4500', clipping=None, suffix='_sm_detrend')
    fig.suptitle(r"SpecMatch-Emp Results, $T_{eff} < 4500$ K, Method: Linear Combination, Detrended"+"\n"+title, fontsize=18)
    plt.tight_layout(rect=[0,0.03,1,0.95])
    plt.savefig(os.path.join(outdir, prefix+"_diagnostic_cool_detrend.png"))
    plt.close(fig)

    fig = plt.figure(figsize=(15,12))
    plots.diagnostic_plots(lib, query='Teff > 4500', clipping=None, suffix='_sm_detrend')
    fig.suptitle(r"SpecMatch-Emp Results, $T_{eff} > 4500$ K, Method: Linear Combination, Detrended"+"\n"+title, fontsize=18)
    plt.tight_layout(rect=[0,0.03,1,0.95])
    plt.savefig(os.path.join(outdir, prefix+"_diagnostic_hot_detrend.png"))
    plt.close(fig)



if __name__ == '__main__':
    # Argument parser
    psr = ArgumentParser(description="Produce plots for analysis of match results")
    psr.add_argument('library', type=str, help="Path to library h5 file")
    psr.add_argument('results', type=str, help="Format string for path to results csv files")
    psr.add_argument('outdir', type=str, help="Path to output directory")
    psr.add_argument('prefix', type=str, help="String to prefix to output files")
    psr.add_argument('title', type=str, help="String to add to title")
    psr.add_argument('wavelengths', type=int, nargs='*', help="Wavelength regions to average")
    args = psr.parse_args()

    mpl.rcParams['figure.dpi'] = 200


    if not os.path.isfile(args.library):
        print("Could not find {0}".format(args.library))
        sys.exit()

    main(args.library, args.results, args.outdir, args.prefix, args.title, args.wavelengths)