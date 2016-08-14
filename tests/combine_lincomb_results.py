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

    min_num = 2
    max_num = 9
    for name in names:
        # combine results as one row for each star.
        res_star = pd.DataFrame()
        for num_best in range(min_num, max_num):
            res_path = os.path.join(args.resdir, '{0}/{0}_lincomb{1:d}.csv'.format(name, num_best))

            # individual dataframe for each and lincomb number
            res_star_num = pd.read_csv(res_path, index_col=0)

            if res_star.empty:
                res_star = res_star_num
                res_star['targ_idx'] = res_star.iloc[0].name
            else:
                # join on columns
                cols = [l for l in res_star_num.columns if 'lincomb{0:d}'.format(num_best) in l]
                res_star[cols] = res_star_num[cols]

        # save result for each star
        outpath_star = os.path.join(args.resdir, '/{0}/{0}_lincomb.csv'.format(name))
        res_star.to_csv(outpath_star)

        # concatenate each star to global dataframe
        if res_global.empty:
            res_global = res_star
        else:
            res_global = pd.concat((res_global, res_star), ignore_index=True)

    # save global results
    res_global.reset_index(inplace=True)
    outpath_global = os.path.join(args.resdir, 'lincomb_results.csv')
    res_global.to_csv(outpath_global)
