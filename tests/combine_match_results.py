#!/usr/bin/env python

import pandas as pd
import os
import sys
import re
from argparse import ArgumentParser

from specmatchemp import library

if __name__ == '__main__':
    # Argument parser
    psr = ArgumentParser(description="Combine match results")
    psr.add_argument('libpath', type=str, help="Path to library file")
    psr.add_argument('resdir', type=str, help="Path to results directory")
    psr.add_argument('-s', '--suffix', type=str, default="", help="Suffix to append to results files")
    args = psr.parse_args()

    if not os.path.isdir(args.resdir):
        print("Could not find folder at {0}".format(args.resdir))
        sys.exit(1)

    lib = library.read_hdf(args.libpath, wavlim='none')

    # get names of stars
    names = list(lib.library_params.cps_name)

    # global dataframe to store 
    res_global = pd.DataFrame()

    for name in names:
        res_path = os.path.join(args.resdir, '{0}/{0}_match.csv'.format(name))

        # individual dataframe for each star
        res_star = pd.read_csv(res_path, index_col=0)

        # concatenate rows to global dataframe
        targ_idx = lib.library_params[lib.library_params.cps_name.str.contains('^'+name+'$')].index[0]
        res_star['targ_idx'] = targ_idx
        res_star['ref_idx'] = res_star.index

        if res_global.empty:
            res_global = res_star
        else:
            res_global = pd.concat((res_global, res_star), ignore_index=True)

    # save global results
    res_global.reset_index(inplace=True)
    outpath_global = os.path.join(args.resdir, 'match_results.csv')
    res_global.to_csv(outpath_global)
