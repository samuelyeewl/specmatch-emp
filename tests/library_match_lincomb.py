#!/usr/bin/env python

import numpy as np
import pandas as pd 
from specmatchemp import match
from specmatchemp import library
import lmfit
import itertools

import os
import sys
from argparse import ArgumentParser

import time

if __name__ == '__main__':
    # Argument parser
    psr = ArgumentParser(description="Cross-match the SpecMatch-Emp library with itself")
    psr.add_argument('library', type=str, help="Path to library")
    psr.add_argument('match_results', type=str, help="Path to match results")
    psr.add_argument('outpath', type=str, help="Path to output file")
    psr.add_argument('num_best', type=int, default=5, help="Number of best matches to consider")
    psr.add_argument('min_w', type=float, help="Minimum wavelength")
    psr.add_argument('length_w', type=float, help="Length of wavelength region")
    args = psr.parse_args()

    # check if library file exists
    if not os.path.isfile(args.library):
        print("Could not find file at {0}".format(args.library))
        sys.exit()
    # read in library
    lib = library.read_hdf(args.library, wavlim=(args.min_w, args.min_w+args.length_w))

    # check if match results file exists
    if not os.path.isfile(args.match_results):
        print("Could not find file at {0}".format(args.match_results))
        sys.exit()
    # read in matches
    matches = pd.DataFrame.from_csv(args.match_results, index_col=0)
    grouped_matches = matches.groupby('targ_idx')

    num_best = args.num_best

    cols = ['targ_idx', 'ref_idxs', 'coeffs', 'chi_squared', 'fit_params']
    chisq_results = []

    # progress meter initialization
    total_matches = len(lib)
    count = 0
    prev_perc = -1 

    # for (param, spec) in [lib[0]]:
    for (param, spec) in lib:
        # progress meter
        perc = int(count/total_matches*100)
        if perc > prev_perc:
            print("{0:d}% complete".format(perc))
            prev_perc += 1
        count += 1

        # get best matches
        best_matches = grouped_matches.get_group(param.lib_index).sort_values(by='chi_squared')
        ref_idxs = np.array(best_matches.head(num_best).ref_idx)
        spec_refs = lib.library_spectra[ref_idxs]

        # get rotational broadening
        vsini = []
        for i in range(num_best):
            params = lmfit.Parameters()
            params.loads(best_matches.iloc[i].fit_params)
            vsini.append(params['vsini'].value)
        vsini = np.array(vsini)

        mt = match.MatchLincomb(lib.wav, spec, spec_refs, vsini)
        mt.best_fit()
        
        # save results
        coeffs = []
        # get coefficients
        for i in range(num_best):
            p = 'coeff_{0:d}'.format(i)
            coeffs.append(mt.best_params[p].value)
        coeffs = np.array(coeffs)

        res = [param.lib_index, ref_idxs, coeffs, mt.best_chisq, mt.best_params.dumps()]  
        chisq_results.append(res)

    df = pd.DataFrame(chisq_results, columns=cols)
    df['targ_idx'] = df['targ_idx'].astype(np.int64)
    df.to_csv(args.outpath)
