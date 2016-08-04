"""
@filename plot_match.py

Plot shift results
"""

import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
import lmfit

import re
import os
from argparse import ArgumentParser

from specmatchemp import library
from specmatchemp import analysis
from specmatchemp import match
from specmatchemp.plotting import plots
from specmatchemp.io import specmatchio

def plot_bestmatch_comparison(res, targ_idx, num_best, cscol, distcol='dist'):
    gs = gridspec.GridSpec(2,2)
    plt.subplot(gs[0])
    plots.bestmatch_comparison_plot(res, 'Teff', 'radius', num_best, cscol, distcol)
    ax = plt.gca()
    ax.set_yscale('log')
    plt.ylim((0.1, 16))
    plots.reverse_x()
    plt.xlabel(r'$T_{eff}$ (K)')
    plt.ylabel(r'$R\ (R_\odot)$')

    plt.subplot(gs[1])
    plots.bestmatch_comparison_plot(res, 'Teff', 'radius', num_best, cscol, distcol)
    ax = plt.gca()
    ax.set_yscale('log')
    ax.yaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
    targ_param = res.loc[targ_idx]
    plt.xlim((targ_param.Teff-200, targ_param.Teff+200))
    plt.ylim((targ_param.radius-0.3*targ_param.radius, targ_param.radius+0.3*targ_param.radius))
    plots.reverse_x()
    plt.xlabel(r'$T_{eff}$ (K)')
    plt.ylabel(r'$R\ (R_\odot)$')

    plt.subplot(gs[2])
    plots.bestmatch_comparison_plot(res, 'feh', 'radius', num_best, cscol, distcol)
    ax = plt.gca()
    ax.set_yscale('log')
    plt.ylim((0.1, 16))
    plt.xlabel(r'$[Fe/H]$ (dex)')
    plt.ylabel(r'$R\ (R_\odot)$')
    
    plt.subplot(gs[3])
    plots.bestmatch_comparison_plot(res, 'feh', 'radius', num_best, cscol, distcol)
    ax = plt.gca()
    ax.set_yscale('log')
    ax.yaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
    targ_param = res.loc[targ_idx]
    plt.xlim((targ_param.feh-0.5, targ_param.feh+0.5))
    plt.ylim((targ_param.radius-0.3*targ_param.radius, targ_param.radius+0.3*targ_param.radius))
    plt.xlabel(r'$[Fe/H]$ (dex)')
    plt.ylabel(r'$R\ (R_\odot)$')
    
    fig = plt.gcf()
    gs.tight_layout(fig, rect=[0, 0, 1, 0.95])


def plot_chi_squared(res, targ_idx, suffix):
    """Wrapper function"""
    plots.chi_squared_plot(res, targ_idx, suffix)


def plot_best_matches(res, targ_idx, lib, wl, cscol='chi_squared', fitcol='fit_params', num_best=5):
    targ_spec = lib.library_spectra[targ_idx]
    labels = {}
    labels['targ_label'] = lib.library_params.loc[targ_idx].cps_name

    res.sort_values(by=cscol, inplace=True)

    for i in range(num_best):
        ref_idx = res.iloc[i].ref_idx
        ref_spec = lib.library_spectra[ref_idx]
        labels['ref_label'] = lib.library_params.loc[ref_idx].cps_name
        labels['res_label'] = r'$\chi^2$ = '+'{0:.3f}'.format(res.iloc[i][cscol])

        fit_params = lmfit.Parameters()
        fit_params.loads(res.iloc[i][fitcol])
        labels['mod_label'] = r'$v \sin i$ = '+'{0:.2f}'.format(fit_params['vsini'].value)

        mt = match.Match(lib.wav, targ_spec, ref_spec)
        mt.create_model(fit_params)

        if i == 0:
            ax = plt.subplot(num_best, 1, i+1)
            shared_ax = ax
        else:
            ax = plt.subplot(num_best, 1, i+1, sharex=shared_ax)
        
        ax.set_xlim((wl+50, wl+75))
        ax.set_title('Best Match {0:d}'.format(i+1), fontsize=12)
        
        plots.plot_match(mt, labels=labels)
        
        if i != num_best-1:
            # Hide all but the lowest axis labels
            plt.setp(ax.get_xticklabels(), visible=False)
        else:
            plt.xlabel('Wavelength (Angstroms)')

    plt.tight_layout(rect=[0, 0, 1, 0.95])


def main(cps_name, respath, outpath, libpath, num_best):
    lib = library.read_hdf(libpath)
    res = pd.read_csv(respath, index_col=0)
    targ_idx = lib.get_index(cps_name)
    res.loc[:,'targ_idx'] = targ_idx
    res.loc[:,'ref_idx'] = res.index

    # calculate distance 
    res.loc[:,'dist'] = res.apply(analysis.dist, args=(res.loc[targ_idx],),axis=1)

    # get different wavelength regions
    cscols = [c for c in list(res.columns) if re.search('chi_squared_\d+', c)]
    fitcols = [c for c in list(res.columns) if re.search('fit_params_\d+', c)]
    wls = [re.search('chi_squared_(\d+)$', c).group(1) for c in cscols]
    wls = map(int, wls)


    with PdfPages(outpath) as pdf:
        for (wl, cscol, fitcol) in zip(wls, cscols, fitcols):
            fig = plt.figure(figsize=(12,8))
            plt.suptitle('Chi-squared surfaces for star {0}'.format(cps_name) \
                +'\nWavelength Region: {0:d} - {1:d} A'.format(wl, wl+100), fontsize=16)
            plot_chi_squared(res, targ_idx, '_{0:d}'.format(wl))
            pdf.savefig()
            plt.close()

            fig = plt.figure(figsize=(12,10))
            plt.suptitle('Comparison of best matches for star {0}'.format(cps_name) \
                +'\nWavelength Region: {0:d} - {1:d} A'.format(wl, wl+100), fontsize=16)
            plot_bestmatch_comparison(res, targ_idx, num_best, cscol)
            pdf.savefig()
            plt.close()

            fig = plt.figure(figsize=(10,12))
            plot_best_matches(res, targ_idx, lib, wl, cscol, fitcol, num_best=num_best)
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
