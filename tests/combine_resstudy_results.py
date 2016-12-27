#!/usr/bin/env python
"""
@filename combine_resstudy_results.py

Combine results from resolution degradation study
"""

import pandas as pd
import numpy as np
from collections import defaultdict

from specmatchemp import library
from specmatchemp import specmatch

# Noise selection path
lst = '/home/syee/specmatchemp-working/specmatchemp/tests/noise_selection.csv'
lc = '/home/syee/specmatchemp-working/specmatchemp/results/lincomb_results.csv'
res_path = '/home/syee/specmatchemp-working/specmatchemp/tests/resstudy_result.csv'
resdir = '/home/syee/specmatchemp-working/specmatchemp/results/'
RESOLUTION = [50000, 30000, 20000, 10000, 5000, 2500]

if __name__ == '__main__':
    # Read library
    lib = library.read_hdf()
    # Get list of stars which were tested
    starlist = pd.read_csv(lst)
    # Open lincomb results
    lincomb_res_table = pd.read_csv(lc)

    for idx, row in starlist.iterrows():
        star = row['cps_name']
        # Get original lincomb result
        lc_res = lincomb_res_table.query('cps_name == "' + star + '"').iloc[0]
        suf = '_nodetrend_lincomb'
        lincomb_result = {}
        for p in library.Library.STAR_PROPS:
            p_col = p + suf
            starlist.loc[idx, p_col] = lc_res[p_col]
            lincomb_result[p] = lc_res[p_col]

        # Read results for each resolution
        for res in RESOLUTION:
            result_list = {}
            # Read in results for each realization
            sm_name = resdir + '{0}/{0}/{0}_R={1:d}.hdf'.format(star, res)
            sm = specmatch.SpecMatch.read_hdf(sm_name, lib)
            for p in library.Library.STAR_PROPS:
                result_list[p] = sm.results_nodetrend[p] - lincomb_result[p]

            p_col = p + '_R={0:d}'.format(res)
            starlist.loc[idx, p_col] = result_list[p]

    starlist.to_csv(res_path)
