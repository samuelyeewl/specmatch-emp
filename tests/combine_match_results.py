#!/usr/bin/env python

import pandas as pd
import glob
import os
import sys
import re
from argparse import ArgumentParser

if __name__ == '__main__':
    # Argument parser
    psr = ArgumentParser(description="Combine match results")
    psr.add_argument('resdir', type=str, help="Path to results directory")
    psr.add_argument('-s', '--suffix', type=str, default="", help="Suffix to append to results files")
    args = psr.parse_args()

    if not os.path.isdir(args.resdir):
        print("Could not find folder at {0}".format(args.resdir))
        sys.exit(1)

    # get names of stars
    dirs = glob.glob(os.path.join(args.resdir, '*/'))
    names = [os.path.basename(os.path.normpath(d)) for d in dirs]

    # global dataframe to store 
    res_global = pd.DataFrame()
    for name, sdir in zip(names, dirs):
        # individual dataframe for each star
        res_star = pd.DataFrame()

        # get wavelengths used
        files = glob.glob(os.path.join(sdir, name+'_*_match.csv'))
        files = [f for f in files if re.search(str(name)+r'_\d+_match.csv', f)]
        wls = [re.search(name+'_(.+?)_match.csv', f).group(1) for f in files]

        for wl, f in zip(wls, files):
            if res_star.empty:
                res_star = pd.read_csv(f, index_col=0)
                res_star.rename(columns={'lib_index.1':'lib_index'}, inplace=True)
            else:
                # merge on columns
                df = pd.read_csv(f, index_col=0)
                res_star = res_star.join(df['chi_squared_{0} fit_params_{0}'.format(wl).split()])

        # save each star's result
        outpath_star = os.path.join(sdir, name+'_match.csv')
        res_star.to_csv(outpath_star)

        # concatenate rows to global dataframe
        targ_idx = res_star[res_star.cps_name.str.contains('^'+str(name)+'$')].index[0]
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
