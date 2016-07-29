#!/usr/bin/env python
"""
@filename match_library_spectrum.py

Matches a spectrum against the rest of the library, across
all wavelengths in 100 A segments.
Saves results in a csv file.
"""

import os
import numpy as np
import pandas as pd 
from specmatchemp import match
from specmatchemp import library

from argparse import ArgumentParser

def main(libpath, targ_name, outpath, wavlim):
    lib = library.read_hdf(libpath, wavlim=wavlim)
    params = lib.library_params

    targ_idx = lib.get_index(targ_name)
    param_targ, spec_targ = lib[targ_idx]

    # create columns 
    cs_col = 'chi_squared_{0:d}'.format(wavlim[0])
    params.loc[:,cs_col] = np.nan
    fit_col = 'fit_params_{0:d}'.format(wavlim[0])
    params.loc[:,fit_col] = ""

    for param_ref, spec_ref in lib:
        # don't match a spectrum against itself
        if param_ref.cps_name == param_targ.cps_name:
            continue

        # match
        mt = match.Match(lib.wav, spec_targ, spec_ref, opt='nelder')
        mt.best_fit()

        # store results
        ref_idx = param_ref.lib_index
        params.loc[ref_idx,cs_col] = mt.best_chisq
        params.loc[ref_idx,fit_col] = mt.best_params.dumps()

    params.to_csv(outpath)


if __name__ == '__main__':
    # Argument parser
    psr = ArgumentParser(description="Cross-match the SpecMatch-Emp library with itself")
    psr.add_argument('library', type=str, help="Path to library")
    psr.add_argument('cps_name', type=str, help="Target name")
    psr.add_argument('outdir', type=str, help="Path to output directory")
    psr.add_argument('minw', type=int, help="Minimum wavelength (A)")
    psr.add_argument('lengthw', type=int, help="Length of wavlength region (A)")
    psr.add_argument('-s', '--suffix', type=str, default="", help="Suffix to append to results files")
    args = psr.parse_args()

    outdir = os.path.join(args.outdir, args.cps_name)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = os.path.join(outdir, args.cps_name+'_{0:d}_match{1}.csv'.format(args.minw, args.suffix))
    wavlim = (args.minw, args.minw+args.lengthw)

    main(args.library, args.cps_name, outpath, wavlim)