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
from specmatchemp import specmatch
from specmatchemp import library

from argparse import ArgumentParser

WAVLIM=(5000,6400)
WAVSTEP=100

def main(libpath, targ_name, outpath):
    lib = library.read_hdf(libpath)
    targ_idx = lib.get_index(targ_name)
    targ_param, targ_spec = lib[targ_idx]

    results = lib.library_params.copy()
    
    for wl in range(WAVLIM[0], WAVLIM[1], WAVSTEP):
        print("Matching region {0:d} - {1:d}".format(wl, wl+100))
        sm = specmatch.SpecMatch(targ_spec, lib, (wl,wl+WAVSTEP))
        sm.match(lincomb=False, ignore=targ_idx)

        cs_col = 'chi_squared_{0:d}'.format(wl)
        results[cs_col] = sm.match_results['chi_squared']
        fit_col = 'fit_params_{0:d}'.format(wl)
        results[fit_col] = sm.match_results['fit_params']

    results.to_csv(outpath)


if __name__ == '__main__':
    # Argument parser
    psr = ArgumentParser(description="Cross-match the SpecMatch-Emp library with itself")
    psr.add_argument('library', type=str, help="Path to library")
    psr.add_argument('cps_name', type=str, help="Target name")
    psr.add_argument('outdir', type=str, help="Path to output directory")
    psr.add_argument('-s', '--suffix', type=str, default="", help="Suffix to append to results files")
    args = psr.parse_args()

    outdir = os.path.join(args.outdir, args.cps_name)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = os.path.join(outdir, args.cps_name+'_match{0}.csv'.format(args.suffix))

    main(args.library, args.cps_name, outpath)