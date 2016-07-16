#!/usr/bin/env python

import pandas as pd 
from specmatchemp import match
from specmatchemp import library

import os
import sys
from argparse import ArgumentParser

import time

def main(libpath, outpath, targ_idx, min_w, length_w):
    # read in library
    lib = library.read_hdf(libpath, wavlim=(min_w, min_w+length_w))

    cols = ['targ_idx', 'ref_idx', 'chi_squared', 'fit_params']
    chisq_results = []

    param, spec = lib[targ_idx]

    for param_ref, spec_ref in lib:
        # don't match a spectrum against itself
        if param.cps_name == param_ref.cps_name:
            continue
        else:
            mt = match.Match(lib.wav, spec, spec_ref, opt='lm')
            mt.best_fit()

            # record results
            res = [param.lib_index, param_ref.lib_index, mt.best_chisq, mt.best_params.dumps()]
            chisq_results.append(res)

    df = pd.DataFrame(chisq_results, columns=cols)
    df.to_csv(args.outpath)


if __name__ == '__main__':
    # Argument parser
    psr = ArgumentParser(description="Cross-match the SpecMatch-Emp library with itself")
    psr.add_argument('library', type=str, help="Path to library")
    psr.add_argument('outpath', type=str, help="Path to output file")
    psr.add_argument('targ_idx', type=int, help="Target index")
    psr.add_argument('min_w', type=float, help="Minimum wavelength")
    psr.add_argument('length_w', type=float, help="Length of wavelength region")
    args = psr.parse_args()

    # check if library file exists
    if not os.path.isfile(args.library):
        print("Could not find file at {0}".format(args.library))
        sys.exit()

    main(args.library, args.outpath, args.targ_idx, args.min_w, args.length_w)
