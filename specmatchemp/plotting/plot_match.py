#!/usr/bin/env python
"""
@filename plot_match.py

Plot shift results
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

import re
import os
from argparse import ArgumentParser

from specmatchemp import library
from specmatchemp import analysis
from specmatchemp.specmatch import SpecMatch


def main(cps_name, respath, outpath, libpath, num_best):
    lib = library.read_hdf(libpath)
    res = pd.read_csv(respath, index_col=0)
    targ_idx = lib.get_index(cps_name)
    targ_param, targ_spec = lib.pop(targ_idx)
    res.loc[:,'targ_idx'] = targ_idx
    res.loc[:,'ref_idx'] = res.index

    # calculate distance 
    res.loc[:,'dist'] = res.apply(analysis.dist, args=(targ_param,),axis=1)
    res.sort_values(by='dist', inplace=True)

    # get different wavelength regions
    cscols = [c for c in list(res.columns) if re.search('chi_squared_\d+', c)]
    fitcols = [c for c in list(res.columns) if re.search('fit_params_\d+', c)]
    wls = [re.search('chi_squared_(\d+)$', c).group(1) for c in cscols]
    wls = map(int, wls)

    with PdfPages(outpath) as pdf:
        for (wl, cscol, fitcol) in zip(wls, cscols, fitcols):
            sm = SpecMatch(targ_spec, lib, (wl,wl+100), num_best)
            sm.match_results = res.rename(columns={cscol:'chi_squared', fitcol:'fit_params'})
            sm.match_results.sort_values(by='chi_squared', inplace=True)

            ### chi-squared surface
            fig = plt.figure(figsize=(12,8))
            plt.suptitle('Chi-squared surfaces for star {0}'.format(cps_name) \
                +'\nWavelength Region: {0:d} - {1:d} A'.format(wl, wl+100), fontsize=16)
            sm.plot_chi_squared_surface()
            # plot actual parameters
            axes = fig.axes
            axes[0].axvline(targ_param['Teff'], color='k')
            axes[1].axvline(targ_param['radius'], color='k')
            axes[2].axvline(targ_param['feh'], color='k')
            pdf.savefig()
            plt.close()

            ### best match comparison
            fig = plt.figure(figsize=(12,10))
            plt.suptitle('Comparison of best matches for star {0}'.format(cps_name) \
                +'\nWavelength Region: {0:d} - {1:d} A'.format(wl, wl+100), fontsize=16)
            sm.plot_references(verbose=False)
            axes = fig.axes
            # plot target
            axes[0].plot(targ_param['Teff'], targ_param['radius'], '*', ms=15, color='red', label='Target')
            axes[1].plot(targ_param['Teff'], targ_param['radius'], '*', ms=15, color='red')
            axes[2].plot(targ_param['feh'], targ_param['radius'], '*', ms=15, color='red')
            axes[3].plot(targ_param['feh'], targ_param['radius'], '*', ms=15, color='red')

            # plot actual closest stars
            axes[0].plot(res.head(num_best)['Teff'], res.head(num_best)['radius'], 's', color='magenta', label='Closest Stars')
            axes[1].plot(res.head(num_best)['Teff'], res.head(num_best)['radius'], 's', color='magenta')
            axes[2].plot(res.head(num_best)['feh'], res.head(num_best)['radius'], 's', color='magenta')
            axes[3].plot(res.head(num_best)['feh'], res.head(num_best)['radius'], 's', color='magenta')

            axes[0].legend(numpoints=1, fontsize='small', loc='best')
            pdf.savefig()
            plt.close()

            ### best matches
            fig = plt.figure(figsize=(10,12))
            sm.plot_best_matches(num_best)
            pdf.savefig()
            plt.close()


if __name__ == '__main__':
    psr = ArgumentParser(description="Build the SpecMatch-Emp library from the various catalogs")
    psr.add_argument('name', type=str, help="CPS name of star")
    psr.add_argument('resdir', type=str, help="Directory of results")
    psr.add_argument('libpath', type=str, help="Path to library")
    psr.add_argument('-n', '--num_best', type=int, default=5, help="Number of best matches to plot")
    psr.add_argument('-s', '--suffix', type=str, default="", help="Suffix on results file")
    args = psr.parse_args()

    respath = os.path.join(args.resdir+args.name+"/"+args.name+args.suffix+"_match.csv")
    outpath = os.path.join(args.resdir+args.name+"/"+args.name+args.suffix+"_match.pdf")

    main(args.name, respath, outpath, args.libpath, args.num_best)
