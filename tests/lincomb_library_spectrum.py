#!/usr/bin/env python
"""
@filename lincomb_library_spectrum.py

Uses the linear combination approach to find the best match for a library spectrum
"""

import numpy as np
import pandas as pd 
from specmatchemp import specmatch
from specmatchemp import match
from specmatchemp import library
import lmfit
import itertools
import json

import os
import sys
from argparse import ArgumentParser

WAVLIM=(5000,6400)
WAVSTEP=100

def main(libpath, targ_name, respath, outpath, num_best):
    lib = library.read_hdf(libpath)
    targ_idx = lib.get_index(targ_name)
    targ_param, targ_spec = lib[targ_idx]

    res_match = pd.read_csv(respath, index_col=0)
    res_lincomb = targ_param.to_frame().transpose()

    for wl in range(WAVLIM[0], WAVLIM[1], WAVSTEP):
        matchcs_col = 'chi_squared_{0:d}'.format(wl)
        res_match.sort_values(by=matchcs_col, inplace=True)
        matchfit_col = 'fit_params_{0:d}'.format(wl)

        sm = specmatch.SpecMatch(targ_spec, lib, (wl,wl+100), num_best)
        sm.match(match_results=res_match.rename(\
            columns={matchcs_col:'chi_squared', matchfit_col:'fit_params'}))

        bestcs_col = 'best_cs_lincomb{0:d}_{1:d}'.format(num_best, wl)
        bestcs = np.array(res_match.head(num_best)[matchcs_col])
        res_lincomb[bestcs_col] = json.dumps(bestcs.tolist())

        refidxs_col = 'ref_idxs_lincomb{0:d}_{1:d}'.format(num_best, wl)
        ref_idxs = np.array(res_match.head(num_best).index)
        res_lincomb[refidxs_col] = json.dumps(ref_idxs.tolist())

        coeffs_col = 'coeffs_lincomb{0:d}_{1:d}'.format(num_best, wl)
        coeffs = np.array(match.get_lincomb_coeffs(sm.mt_lincomb.best_params))
        res_lincomb[coeffs_col] = json.dumps(coeffs.tolist())

        cs_col = 'chi_squared_lincomb{0:d}_{1:d}'.format(num_best, wl)
        res_lincomb[cs_col] = sm.mt_lincomb.best_chisq

        fit_col = 'fit_params_lincomb{0:d}_{1:d}'.format(num_best, wl)
        res_lincomb[fit_col] = sm.mt_lincomb.best_params.dumps()

    res_lincomb.to_csv(outpath)


if __name__ == '__main__':
    psr = ArgumentParser(description="Cross-match the SpecMatch-Emp library with itself")
    psr.add_argument('library', type=str, help="Path to library")
    psr.add_argument('cps_name', type=str, help="Target name")
    psr.add_argument('resdir', type=str, help="Path to results directory")
    psr.add_argument('num_best', type=int, help="Number of best matches to consider")
    psr.add_argument('-s', '--suffix', type=str, default="", help="Suffix to append to results files")
    args = psr.parse_args()

    if not os.path.isfile(args.library):
        print("Could not find library at {0}".format(args.library))
        sys.exit(1)

    outdir = os.path.join(args.resdir, args.cps_name)
    if not os.path.exists(outdir):
        print("Star {0} was not found".format(args.cps_name))
        sys.exit(1)

    respath = os.path.join(outdir, args.cps_name \
        + '_match{0}.csv'.format(args.suffix))
    outpath = os.path.join(outdir, args.cps_name \
        +'_lincomb{0:d}{1}.csv'.format(args.num_best, args.suffix))

    main(args.library, args.cps_name, respath, outpath, args.num_best)

