#!/usr/bin/env python
"""
@filename lincomb_library_spectrum.py

Uses the linear combination approach to find the best match for a library spectrum
"""

import numpy as np
import pandas as pd 
from specmatchemp import match
from specmatchemp import library
import lmfit
import itertools
import json

import os
import sys
from argparse import ArgumentParser

def main(libpath, targ_name, respath, outpath, num_best, wavlim):
    lib = library.read_hdf(libpath, wavlim=wavlim)
    targ_idx = lib.get_index(targ_name)
    param_targ, spec_targ = lib[targ_idx]

    res_match = pd.read_csv(respath, index_col=0)
    res_match.sort_values(by='chi_squared_{0}'.format(wavlim[0]), inplace=True)

    res_lincomb = param_targ.to_frame().transpose()
    bestcs_col = 'best_cs_lincomb{0:d}_{1:d}'.format(num_best, wavlim[0])
    refidxs_col = 'ref_idxs_lincomb{0:d}_{1:d}'.format(num_best, wavlim[0])
    coeffs_col = 'coeffs_lincomb{0:d}_{1:d}'.format(num_best, wavlim[0])
    cs_col = 'chi_squared_lincomb{0:d}_{1:d}'.format(num_best, wavlim[0])
    fit_col = 'fit_params_lincomb{0:d}_{1:d}'.format(num_best, wavlim[0])

    # get best matches
    bestcs = np.array(res_match.head(num_best)['chi_squared_{0:d}'.format(wavlim[0])])
    res_lincomb[bestcs_col] = json.dumps(bestcs.tolist())
    ref_idxs = np.array(res_match.head(num_best).index)
    res_lincomb[refidxs_col] = json.dumps(ref_idxs.tolist())
    # get reference spectra
    spec_refs = lib.library_spectra[ref_idxs]
    # get vsini
    vsini = []
    for i in range(num_best):
        params = lmfit.Parameters()
        params.loads(res_match.iloc[i]['fit_params_{0:d}'.format(wavlim[0])])
        vsini.append(params['vsini'].value)
    vsini = np.array(vsini)

    # perform fit
    mt = match.MatchLincomb(lib.wav, spec_targ, spec_refs, vsini)
    mt.best_fit()

    # save result
    coeffs = np.array(match.get_lincomb_coeffs(mt.best_params))
    res_lincomb[coeffs_col] = json.dumps(coeffs.tolist())
    res_lincomb[cs_col] = mt.best_chisq
    res_lincomb[fit_col] = mt.best_params.dumps()

    res_lincomb.to_csv(outpath)


if __name__ == '__main__':
    psr = ArgumentParser(description="Cross-match the SpecMatch-Emp library with itself")
    psr.add_argument('library', type=str, help="Path to library")
    psr.add_argument('cps_name', type=str, help="Target name")
    psr.add_argument('resdir', type=str, help="Path to results directory")
    psr.add_argument('num_best', type=int, help="Number of best matches to consider")
    psr.add_argument('minw', type=int, help="Minimum wavelength (A)")
    psr.add_argument('lengthw', type=int, help="Length of wavlength region (A)")
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
        + '_{0:d}_match{1}.csv'.format(args.minw, args.suffix))
    outpath = os.path.join(outdir, args.cps_name \
        +'_lincomb{0:d}_{1}_match{2}.csv'.format(args.num_best, args.minw, args.suffix))
    wavlim = (args.minw, args.minw+args.lengthw)

    main(args.library, args.cps_name, respath, outpath, args.num_best, wavlim)

