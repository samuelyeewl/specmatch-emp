#!/usr/bin/env python

import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from specmatch_io import *
from specmatchemp import match
import lmfit

if __name__ == '__main__':
    lib = pd.read_csv('~/SpecMatch-Emp/starswithspectra.csv',index_col=0)
    lib = lib.convert_objects(convert_numeric=True)
    lib["chi-squared"] = np.nan
    lib["params"] = ''

    specref_path = '/Users/samuel/Dropbox/SpecMatch-Emp/nso/nso_std.fits'
    specref, headerref = read_as_dataframe(specref_path)
    query = '6000 < w < 6100'
    specref = specref.query(query)

    for index, row in lib.iterrows():
        spec_path = '/Users/samuel/SpecMatch-Emp/lib/'+row['obs']+'_adj.fits'
        try:
            spec, header = read_as_dataframe(spec_path)
        except Exception as e:
            continue
        mt = match.Match(spec, specref)
        mt.best_fit()
        lib.loc[index, "chi-squared"] = mt.best_chisq
        print(mt.best_chisq)
        lib.loc[index, "params"] = mt.best_params.dumps()

    lib.to_csv('./chi-squared.csv')