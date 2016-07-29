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


import os
from argparse import ArgumentParser

from specmatchemp import library
from specmatchemp import analysis
from specmatchemp.plotting import plots
from specmatchemp.io import specmatchio

def plot_bestmatch_comparison(targ_idx, res, num_best, cscol, distcol='dist'):
    gs = gridspec.GridSpec(2,2)
    plt.subplot(gs[0])
    plots.bestmatch_comparison_plot(res, 'Teff', 'logg', num_best, cscol, distcol)
    plots.reverse_x()
    plots.reverse_y()
    plt.subplot(gs[1])
    plots.bestmatch_comparison_plot(res, 'Teff', 'logg', num_best, cscol, distcol)
    targ_param = res.loc[targ_idx]
    plt.xlim((targ_param.Teff-200, targ_param.Teff+200))
    plt.ylim((targ_param.logg-0.5, targ_param.logg+0.5))
    plots.reverse_x()
    plots.reverse_y()

    plt.subplot(gs[2])
    plots.bestmatch_comparison_plot(res, 'feh', 'logg', num_best, cscol, distcol)
    plots.reverse_y()
    plt.subplot(gs[3])
    plots.bestmatch_comparison_plot(res, 'feh', 'logg', num_best, cscol, distcol)
    targ_param = res.loc[targ_idx]
    plt.xlim((targ_param.feh-0.5, targ_param.feh+0.5))
    plt.ylim((targ_param.logg-0.5, targ_param.logg+0.5))
    plots.reverse_y()
    fig = plt.gcf()
    gs.tight_layout(fig, rect=[0, 0, 1, 0.95])


def main(cps_name, respath, outpath, libpath, num_best, wavlim):
    lib = library.read_hdf(libpath, wavlim=wavlim)
    res = pd.read_csv(respath, index_col=0)
    targ_idx = lib.get_index(cps_name)
    res.loc[:,'targ_idx'] = targ_idx
    res.loc[:,'ref_idx'] = res.index

    # calculate distance 
    res.loc[:,'dist'] = res.apply(analysis.dist, args=(res.loc[targ_idx],),axis=1)

    with PdfPages(outpath) as pdf:
        fig = plt.figure(figsize=(12,8))
        plt.suptitle('Chi-squared surfaces for star {0}'.format(cps_name), fontsize=16)
        plots.chi_squared_plot(res, targ_idx, '_{0:d}'.format(wavlim[0]))
        pdf.savefig()
        plt.close()

        fig = plt.figure(figsize=(12,10))
        plt.suptitle('Comparison of best matches for star {0}'.format(cps_name), fontsize=16)
        plot_bestmatch_comparison(targ_idx, res, num_best, 'chi_squared_{0:d}'.format(wavlim[0]))
        pdf.savefig()
        plt.close()



if __name__ == '__main__':
    psr = ArgumentParser(description="Build the SpecMatch-Emp library from the various catalogs")
    psr.add_argument('name', type=str, help="CPS name of star")
    psr.add_argument('resdir', type=str, help="Directory of results")
    psr.add_argument('libpath', type=str, help="Path to library")
    psr.add_argument('num_best', type=int, help="Number of best matches to plot")
    psr.add_argument('minw', type=int, help="Minimum wavelength (A)")
    psr.add_argument('lengthw', type=int, help="Length of wavlength region (A)")
    psr.add_argument('-s', '--suffix', type=str, default="", help="Suffix on results file")
    args = psr.parse_args()

    respath = os.path.join(args.resdir+args.name+"/"+args.name+args.suffix+"_match.csv")
    outpath = os.path.join(args.resdir+args.name+"/"+args.name+args.suffix+"_match.pdf")
    wavlim = (args.minw, args.minw+args.lengthw)
    main(args.name, respath, outpath, args.libpath, args.num_best, wavlim)