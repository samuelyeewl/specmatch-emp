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

def generate_sm_values(lib, results, method='lincomb', suffix='_sm'):
    params = lib.library_params

    # Use closest matched spectra as sm values
    if method == 'lincomb':
        results.set_index('targ_idx', inplace=True)
        for p in library.STAR_PROPS:
            psm = p+suffix
            params[psm] = params.lib_index.apply(lambda i: \
                lincomb_props(lib.library_params, p, results.loc[i].ref_idxs, results.loc[i].coeffs))
    return params

def lincomb_props(library_params, prop, idxs, coeffs):
    sm_prop = 0
    for i in range(len(idxs)):
        lib_prop = library_params.loc[idxs[i], prop]
        sm_prop += lib_prop*coeffs[i]
    return sm_prop

def main(library_path, results_path, outdir, prefix, title):
    lib = library.read_hdf(library_path)
    results = pd.DataFrame.from_csv(results_path, index_col=0)

    results['ref_idxs'] = results.ref_idxs.apply(lambda s: \
        np.array(s.strip('[]').split()).astype(int))
    results['coeffs'] = results.coeffs.apply(lambda s: \
        np.array(s.strip('[]').split()).astype(float))

    lib.library_params = generate_sm_values(lib, results, method='lincomb', suffix='_sm')

    fig = plt.figure(figsize=(15,12))
    plots.diagnostic_plots(lib, clipping=None, suffix='_sm')
    fig.suptitle("SpecMatch-Emp Results, Method: Linear Combination\n"+title, fontsize=18)
    plt.tight_layout(rect=[0,0.03,1,0.95])
    plt.savefig(os.path.join(outdir, prefix+"_diagnostic_lincomb.png"))
    plt.close(fig)

    fig = plt.figure(figsize=(15,12))
    plots.diagnostic_plots(lib, query='Teff < 4500', clipping=None, suffix='_sm')
    fig.suptitle(r"SpecMatch-Emp Results, $T_{eff} < 4500$ K, Method: Best Match"+"\n"+title, fontsize=18)
    plt.tight_layout(rect=[0,0.03,1,0.95])
    plt.savefig(os.path.join(outdir, prefix+"_diagnostic_cool_lincomb.png"))
    plt.close(fig)



if __name__ == '__main__':
    # Argument parser
    psr = ArgumentParser(description="Produce plots for analysis of match results")
    psr.add_argument('library', type=str, help="Path to library h5 file")
    psr.add_argument('results', type=str, help="Path to results csv file")
    psr.add_argument('outdir', type=str, help="Path to output directory")
    psr.add_argument('prefix', type=str, help="String to prefix to output files")
    psr.add_argument('title', type=str, help="String to add to title")
    args = psr.parse_args()

    mpl.rcParams['figure.dpi'] = 200


    if not os.path.isfile(args.library):
        print("Could not find {0}".format(args.library))
        sys.exit()
    if not os.path.isfile(args.results):
        print("Could not find {0}".format(args.results))
        sys.exit()

    main(args.library, args.results, args.outdir, args.prefix, args.title)