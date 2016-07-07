#!/usr/bin/env python

import pandas as pd 
from specmatchemp import match
from specmatchemp import library
import lmfit
import itertools

WAVLIM = (6000, 6100)

if __name__ == '__main__':
    # read in library
    lib = library.read_hdf('/Users/samuel/SpecMatch-Emp/lib/library.h5',wavlim=WAVLIM)

    cols = ['targ_idx', 'ref_idx', 'chi_squared', 'fit_params']
    chisq_results = []
    # total
    # total_matches = len(lib)**2
    total_matches = len(lib)
    # counter
    i = 0
    prev_perc = -1 
    # for (param, spec), (param_ref, spec_ref) in itertools.product(lib, repeat=2):
    param, spec = lib[0]
    for (param_ref, spec_ref) in lib:
        # progress meter
        # perc = int(i/total_matches*100)
        # if perc > prev_perc:
        #     print("{0:d}% complete".format(perc))
        #     prev_perc = perc
        # i+=1

        # don't match a spectrum against itself
        if param.cps_name == param_ref.cps_name:
            continue
        else:
            mt = match.Match(lib.wav, spec[0], spec[1], spec_ref[0], spec_ref[1])
            mt.best_fit()

            res = [param.lib_index, param_ref.lib_index, \
                    mt.best_chisq, mt.best_params.dumps()]
            chisq_results.append(res)

    df = pd.DataFrame(chisq_results, columns=cols)
    df.to_csv('/Users/samuel/SpecMatch-Emp/tests/chi_squared_results_noerr.csv')
