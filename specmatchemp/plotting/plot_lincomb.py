#!/usr/bin/env python
"""
@filename plot_lincomb.py

Plot lincomb results
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import lmfit

import re
import os
import json
from argparse import ArgumentParser

from specmatchemp import library
from specmatchemp import analysis
from specmatchemp import match
from specmatchemp import plots
from specmatchemp.specmatch import SpecMatch

AVG_PROPS = ['Teff', 'feh', 'radius'] 
WL_AVG = [5000, 5800]

def main(cps_name, respath, matchpath, outpath, libpath, num_best):
    lib = library.read_hdf(libpath)
    res = pd.read_csv(respath, index_col=0)
    targ_idx = lib.get_index(cps_name)
    targ_param, targ_spec = lib[targ_idx]
    res.loc[:,'targ_idx'] = targ_idx

    res_match = pd.read_csv(matchpath, index_col=0)

    cscols = [c for c in list(res.columns) if re.search('chi_squared_lincomb{0:d}_\d+'.format(num_best), c)]
    wls = list(set([re.search('chi_squared_lincomb{0:d}_(\d+)$'.format(num_best), c).group(1) for c in cscols]))
    wls = list(map(int, wls))
    wls.sort()

    wls_avg = [wl for wl in wls if wl > WL_AVG[0] and wl < WL_AVG[1]]
    num_wls = len(wls_avg)

    with PdfPages(outpath) as pdf:
        for wl in wls:
            match_cscol = 'chi_squared_{0:d}'.format(wl)
            match_fitcol = 'fit_params_{0:d}'.format(wl)
            sm = SpecMatch(targ_spec, lib, (wl,wl+100), num_best)
            sm.match_results = res_match.rename(columns={match_cscol:'chi_squared', match_fitcol:'fit_params'})
            sm.match_results.sort_values(by='chi_squared', inplace=True)

            bestcs_col = 'best_cs_lincomb{0:d}_{1:d}'.format(num_best, wl)
            bestcs = json.loads(res.iloc[0][bestcs_col])

            refidxs_col = 'ref_idxs_lincomb{0:d}_{1:d}'.format(num_best, wl)
            refidxs = json.loads(res.iloc[0][refidxs_col])
            sm.ref_idxs = refidxs
            refspecs = lib.get_spectrum(refidxs)
            refspecs = [r.cut(wl,wl+100) for r in refspecs]

            coeffs_col = 'coeffs_lincomb{0:d}_{1:d}'.format(num_best, wl)
            coeffs = json.loads(res.iloc[0][coeffs_col])
            sm.coeffs = coeffs

            cs_col = 'chi_squared_lincomb{0:d}_{1:d}'.format(num_best, wl)
            cs = res.iloc[0][cs_col]

            fit_col = 'fit_params_lincomb{0:d}_{1:d}'.format(num_best, wl)
            fit_params = lmfit.Parameters()
            fit_params.loads(res.iloc[0][fit_col])
            vsini = match.get_vsini(fit_params)

            mt = match.MatchLincomb(targ_spec.cut(wl,wl+100), refspecs, vsini)
            mt.load_params(fit_params)
            mt.best_chisq = cs
            mt.ref_chisq = bestcs

            sm.mt_lincomb = mt
            sm.get_derived_values()

            ## Plot of library best matches
            fig = plt.figure(figsize=(12,10))
            plt.suptitle('Library spectra used to match star {0}'.format(cps_name) \
                +'\nWavelength Region: {0:d} - {1:d} A'.format(wl, wl+100), fontsize=16)
            sm.plot_references()
            axes = fig.axes
            # plot target
            axes[0].plot(targ_param['Teff'], targ_param['radius'], '*', ms=15, color='red', label='Target')
            axes[1].plot(targ_param['Teff'], targ_param['radius'], '*', ms=15, color='red')
            plt.sca(axes[1])
            xlim = plt.xlim()
            ylim = plt.ylim()
            plots.set_tight_lims([xlim[0], xlim[1], targ_param['Teff']], [ylim[0], ylim[1], targ_param['radius']])
            axes[2].plot(targ_param['feh'], targ_param['radius'], '*', ms=15, color='red')
            axes[3].plot(targ_param['feh'], targ_param['radius'], '*', ms=15, color='red')
            plt.sca(axes[3])
            xlim = plt.xlim()
            ylim = plt.ylim()
            plots.set_tight_lims([xlim[0], xlim[1], targ_param['feh']], [ylim[0], ylim[1], targ_param['radius']])

            axes[0].legend(numpoints=1, fontsize='small', loc='best')
            pdf.savefig()
            plt.close()

            ## PLot of match lincomb
            fig = plt.figure(figsize=(12,6))
            sm.plot_lincomb()
            pdf.savefig()
            plt.close()



if __name__ == '__main__':
    psr = ArgumentParser(description="Build the SpecMatch-Emp library from the various catalogs")
    psr.add_argument('libpath', type=str, help="Path to library")
    psr.add_argument('name', type=str, help="CPS name of star")
    psr.add_argument('resdir', type=str, help="Directory of results")
    psr.add_argument('-n', '--num_best', type=int, default=5, help="Number of best matches to plot")
    psr.add_argument('-s', '--suffix', type=str, default="", help="Suffix on results file")
    args = psr.parse_args()

    respath = os.path.join(args.resdir+args.name+"/"+args.name+args.suffix+"_lincomb{0:d}.csv".format(args.num_best))
    matchpath = os.path.join(args.resdir+args.name+"/"+args.name+args.suffix+"_match.csv".format(args.num_best))
    outpath = os.path.join(args.resdir+args.name+"/"+args.name+args.suffix+"_lincomb{0:d}.pdf".format(args.num_best))

    main(args.name, respath, matchpath, outpath, args.libpath, args.num_best)
