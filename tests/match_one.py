#!/usr/bin/env python

import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from specmatchemp.io.io import *
from specmatchemp import match
import lmfit

import time

def main():
    lib = pd.read_csv('~/SpecMatch-Emp/starswithspectra.csv',index_col=0)
    lib = lib.convert_objects(convert_numeric=True)
    lib = lib[lib.Source.str.contains('Mann')]
    lib.reset_index(drop=True, inplace=True)
    lib["chi-squared"] = np.nan
    lib["params"] = ''
    lib["teff_der"] = np.nan
    lib["radius_der"] = np.nan
    lib["mass_der"] = np.nan
    lib["feh_der"] = np.nan
    l = len(lib.index)

    specs = []
    # read in spectra
    for i, row in lib.iterrows():
        spec_path = '/Users/samuel/SpecMatch-Emp/lib_1/'+row['obs']+'_adj.fits'
        try:
            spec, header = read_as_dataframe(spec_path)
        except Exception as e:
            continue
        specs.append(spec)

    starttime = time.time()
    lib_f = lib[lib.obs.str.contains('rj70.782')]
    # for i, row_targ in [(10, lib.iloc[10])]:
    for i, row_targ in lib_f.iterrows():
        libi = lib.copy()
        spec = specs[i]
        for j, row_ref in lib.iterrows():
            if i == j:
                libi.loc[i, "chi-squared"] = np.inf
                continue

            specref = specs[j]

            mt = match.Match(spec, specref)
            mt.best_fit()
            
            libi.loc[j, "chi-squared"] = mt.best_chisq
            libi.loc[j, "params"] = mt.best_params.dumps()

        libi.to_csv('/Users/samuel/SpecMatch-Emp/lib/'+row_targ['obs']+'_cs.csv')
        # find 8 best matches
        num_best = 8
        besti = libi.sort_values(by="chi-squared").head(num_best)
        
        # find a sort of weighted average (testing only)
        besti.loc[:,'weights'] = (1/besti["chi-squared"])/(1/besti["chi-squared"]).sum()
        lib.loc[i, 'teff_der'] = np.average(besti['Teff'], weights=besti['weights'])
        # print(teff_avg, lib.iloc[i].Teff)
        lib.loc[i, 'radius_der'] = np.average(besti['radius'], weights=besti['weights'])
        # print(radius_avg, lib.iloc[i].radius)
        lib.loc[i, 'mass_der'] = np.average(besti['mass'], weights=besti['weights'])
        # print(mass_avg, lib.iloc[i].mass)
        lib.loc[i, 'feh_der'] = np.average(besti['FeH'], weights=besti['weights'])
        # print(feh_avg, lib.iloc[i].FeH)


    lib.to_csv('/Users/samuel/SpecMatch-Emp/lib/updated_chi-squared.csv')
    endtime = time.time()
    print("Took {0:.4f} seconds".format(endtime-starttime))

if __name__ == '__main__':
    main()