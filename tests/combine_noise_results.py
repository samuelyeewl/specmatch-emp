#!/usr/bin/env python
"""
@filename combine_noise_results.py

Combine results from internal cross-validation of library
"""

import pandas as pd
import numpy as np
from collections import defaultdict

from specmatchemp import library
from specmatchemp import specmatch

# Noise selection path
lst = '/home/syee/specmatchemp-working/specmatchemp/tests/noise_selection.csv'
lc = '/home/syee/specmatchemp-working/specmatchemp/results/lincomb_results.csv'
res_path = '/home/syee/specmatchemp-working/specmatchemp/tests/noise_result.csv'
resdir = '/home/syee/specmatchemp-working/specmatchemp/results/'
NUM_SPEC = 20
SNR_LIST = [120, 80, 40, 20, 10]

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

        # Read results for each SNR
        for snr in SNR_LIST:
            result_lists = defaultdict(list)
            # Read in results for each realization
            for i in range(1, NUM_SPEC+1):
                sm_name = resdir + '{0}/{0}/{0}_snr={1:d}_i={2:d}_sm.hdf'\
                    .format(star, snr, i)
                sm = specmatch.SpecMatch.read_hdf(sm_name, lib)
                for p in library.Library.STAR_PROPS:
                    result_lists[p].append(sm.results[p] - lincomb_result[p])

            # Calculate standard deviation of each property
            for p in library.Library.STAR_PROPS:
                p_col = p + '_snr={0:d}'.format(snr)
                starlist.loc[idx, p_col] = np.std(result_lists[p])

    starlist.to_csv(res_path)
