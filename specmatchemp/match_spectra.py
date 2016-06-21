#!/usr/bin/env python
"""
@filename match_spectra.py

Match two spectra
"""

import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from specmatchemp.io.io import *
from specmatchemp import match
import lmfit


if __name__ == '__main__':
    specref_path = '/Users/samuel/Dropbox/SpecMatch-Emp/nso/nso_std.fits'
    spec2_path = './lib/rj118.463_adj.fits'

    specref, headerref = read_as_dataframe(specref_path)
    query = '6000 < w < 6100'
    specref = specref.query(query)

    spec2, header2 = read_as_dataframe(spec2_path)

    mt = match.Match(spec2, specref)
    mt2 = match.Match(spec2,specref)


    mt.best_fit()
    resid = mt.residual(mt.best_params)
    print(mt.best_chisq)
    print(mt.best_params)
    plt.figure(figsize=(12,6))
    plt.plot(mt.spectra['w'], mt.spectra['s_targ'], label='Target')
    plt.plot(mt.spectra['w'], mt.spectra['s_ref'], label='Reference')
    plt.plot(mt.spectra['w'], mt.spectra['s_mod'], label='Modified reference')
    # plt.plot(mt.spectra['w'], resid, label='Residuals')

    # see what happens if we have no vsini
    params = mt.best_params
    params['vsini'].set(value=1.0)
    mt2.create_model(params)
    # plt.plot(mt2.spectra['w'], mt2.spectra['s_mod'], label='Modified reference, no broadening')

    plt.legend(loc='lower right')
    plt.savefig('fit_broadening.png', dpi=300)
    # plt.show()

    # chisquared = calculate_chi_squared(specref, spec2, 6000, 6100)

    # print(chisquared)

    # lib = pd.read_csv('../starswithspectra.csv',index_col=0)
    # lib = lib.convert_objects(convert_numeric=True)
    # lib["chi-squared"] = np.nan

    # for index, row in lib.iterrows():
    #     spec_path = '../lib/'+row['obs']+'_adj.fits'
    #     try:
    #         spec, header = read_as_dataframe(spec_path)
    #     except Exception as e:
    #         continue
    #     lib.loc[index, "chi-squared"] = calculate_chi_squared(specref, spec, 6000, 6100)

    # lib.to_csv('../chi-squared.csv')
