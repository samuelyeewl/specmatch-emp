#!/usr/bin/env python

import numpy as np
import pandas as pd
import glob
import os
import sys
import re
from argparse import ArgumentParser

if __name__ == '__main__':
    # Argument parser
    psr = ArgumentParser(description="Combine lincomb results")
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

        # lincomb num
        min_num = 2
        max_num = 9
        for i in np.arange(min_num, max_num):
            # get wavelengths used
            files = glob.glob(os.path.join(sdir, name+'_lincomb{0:d}_*_match.csv'.format(i)))
            files = [f for f in files if re.search(str(name)+r'_lincomb\d_\d+_match.csv', f)]
            wls = [re.search(name+r'_lincomb\d_(.+?)_match.csv', f).group(1) for f in files]
            wls = map(int, wls)
            for wl, f in zip(wls, files):
                if res_star.empty:
                    res_star = pd.read_csv(f, index_col=0)
                    res_star.rename(columns={'lib_index.1':'lib_index'}, inplace=True)
                    res_star.cps_name = res_star.cps_name.astype(str)
                else:
                    # merge on columns
                    df = pd.read_csv(f, index_col=0)
                    cols = 'best_cs_lincomb{0:d}_{1:d} ref_idxs_lincomb{0:d}_{1:d} coeffs_lincomb{0:d}_{1:d} chi_squared_lincomb{0:d}_{1:d} fit_params_lincomb{0:d}_{1:d}'\
                    .format(i, wl).split()
                    res_star = res_star.join(df[cols])


        # save each star's result
        outpath_star = os.path.join(sdir, name+'_lincomb.csv')
        res_star.to_csv(outpath_star)

        # concatenate rows to global dataframe
        try:
            targ_idx = res_star[res_star.cps_name.str.contains(name)].index[0]
            res_star['targ_idx'] = targ_idx
            if res_global.empty:
                res_global = res_star
            else:
                res_global = pd.concat((res_global, res_star), ignore_index=True)
        except:
            print(name)


    # save global results
    res_global.reset_index(inplace=True)
    outpath_global = os.path.join(args.resdir, 'lincomb_results.csv')
    res_global.to_csv(outpath_global)



