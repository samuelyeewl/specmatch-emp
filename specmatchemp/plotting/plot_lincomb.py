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
import json

from specmatchemp import library
from specmatchemp import analysis
from specmatchemp import match
from specmatchemp.plotting import plots
from specmatchemp.io import specmatchio

AVG_PROPS = ['Teff', 'feh', 'radius'] 
WL_AVG = [5000, 5800]

def plot_lincomb_refs(res, lib, targ_idx, num, wl):
    refidxs_col = 'ref_idxs_lincomb{0:d}_{1:d}'.format(num, wl)
    ref_idxs = json.loads(res.iloc[0][refidxs_col])
    coeffs_col = 'coeffs_lincomb{0:d}_{1:d}'.format(num, wl)
    coeffs = json.loads(res.iloc[0][coeffs_col])
    coeffs = [r"${0:.3f}$".format(c) for c in coeffs]

    gs = gridspec.GridSpec(2,2)
    plt.subplot(gs[0])
    plots.lincomb_refs_plot(lib.library_params, 'Teff', 'radius', targ_idx, ref_idxs)
    ax = plt.gca()
    ax.set_yscale('log')
    plt.ylim((0.1, 16))
    plots.reverse_x()
    plt.xlabel(r'$T_{eff}$ (K)')
    plt.ylabel(r'$R\ (R_\odot)$')

    plt.subplot(gs[1])
    plots.lincomb_refs_plot(lib.library_params, 'Teff', 'radius', targ_idx, ref_idxs, annot=coeffs, zoom=True, legend=False)
    ax = plt.gca()
    # ax.set_yscale('log')
    ax.yaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    plots.reverse_x()
    plt.xlabel(r'$T_{eff}$ (K)')
    plt.ylabel(r'$R\ (R_\odot)$')

    plt.subplot(gs[2])
    plots.lincomb_refs_plot(lib.library_params, 'feh', 'radius', targ_idx, ref_idxs, legend=False)
    ax = plt.gca()
    ax.set_yscale('log')
    plt.ylim((0.1, 16))
    plt.xlabel(r'$[Fe/H]$ (dex)')
    plt.ylabel(r'$R\ (R_\odot)$')

    plt.subplot(gs[3])
    plots.lincomb_refs_plot(lib.library_params, 'feh', 'radius', targ_idx, ref_idxs, annot=coeffs, zoom=True, legend=False)
    ax = plt.gca()
    # ax.set_yscale('log')
    ax.yaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    plt.xlabel(r'$[Fe/H]$ (dex)')
    plt.ylabel(r'$R\ (R_\odot)$')

    plt.tight_layout(rect=[0, 0, 1, 0.95])

def plot_lincomb_specs(res, lib, num, wl):
    targ_idx = res.iloc[0].name
    # Load in values
    bestcs_col = 'best_cs_lincomb{0:d}_{1:d}'.format(num, wl)
    bestcs = json.loads(res.iloc[0][bestcs_col])
    bestcs = ["{0:.2f}".format(c) for c in bestcs]
    refidxs_col = 'ref_idxs_lincomb{0:d}_{1:d}'.format(num, wl)
    ref_idxs = json.loads(res.iloc[0][refidxs_col])
    coeffs_col = 'coeffs_lincomb{0:d}_{1:d}'.format(num, wl)
    coeffs = json.loads(res.iloc[0][coeffs_col])
    coeffs = ["{0:.3f}".format(c) for c in coeffs]
    cs_col = 'chi_squared_lincomb{0:d}_{1:d}'.format(num, wl)
    cs = res.iloc[0][cs_col]
    fit_col = 'fit_params_lincomb{0:d}_{1:d}'.format(num, wl)
    fit_params = lmfit.Parameters()
    fit_params.loads(res.iloc[0][fit_col])
    vsini = match.get_vsini(fit_params)

    wavidx, = np.where(np.logical_and(lib.wav >= wl, lib.wav < wl+100))
    idxmin = wavidx[0]
    idxmax = wavidx[-1]+1

    mt = match.MatchLincomb(lib.wav[idxmin:idxmax], lib.library_spectra[targ_idx,:,idxmin:idxmax], \
        lib.library_spectra[ref_idxs,:,idxmin:idxmax], vsini)
    mt.create_model(fit_params)

    targ_label = str(lib.library_params.loc[targ_idx, 'cps_name'])
    ref_labels=[]
    for i in range(mt.num_refs):
        ref_labels.append(r'{0}, $v \sin i={1:.2f}$, $\chi^2={2}$, $c_{3:d}={4}$'.format(\
            lib.library_params.loc[ref_idxs[i], 'cps_name'], vsini[i], bestcs[i], i+1, coeffs[i]))
    mod_label = r'$\chi^2 = {0:.2f}$'.format(cs)
    plots.plot_matchlincomb(mt, targ_label, ref_labels, mod_label)

    plt.tight_layout(rect=[0, 0, 1, 0.95])

def main(cps_name, respath, outpath, libpath, num):
    lib = library.read_hdf(libpath)
    res = pd.read_csv(respath, index_col=0)
    targ_idx = lib.get_index(cps_name)

    # get wavelength regions and lincomb numbers used
    cscols = [c for c in list(res.columns) if re.search('chi_squared_lincomb{0:d}_\d+'.format(num), c)]
    wls = list(set([re.search('chi_squared_lincomb{0:d}_(\d+)$'.format(num), c).group(1) for c in cscols]))
    wls.sort()
    wls = list(map(int, wls))
    # number of wavelength regions being averaged over
    wls_avg = [wl for wl in wls if wl > WL_AVG[0] and wl < WL_AVG[1]]
    num_wls = len(wls_avg)

    with PdfPages(outpath) as pdf:
        for wl in wls:
            fig = plt.figure(figsize=(12,10))
            plt.suptitle('References used for linear combinations for star {0}'.format(cps_name) \
                +'\nWavelength Region: {0:d} - {1:d} A'.format(wl, wl+100), fontsize=16)
            plot_lincomb_refs(res, lib, targ_idx, num, wl)
            pdf.savefig()
            plt.close()
            
            fig = plt.figure(figsize=(12,6))
            plot_lincomb_specs(res, lib, num, wl)
            plt.xlabel('Wavelength (Angstroms)')
            plt.xlim((wl+50, wl+75))
            pdf.savefig()
            plt.close()


if __name__ == '__main__':
    psr = ArgumentParser(description="Build the SpecMatch-Emp library from the various catalogs")
    psr.add_argument('name', type=str, help="CPS name of star")
    psr.add_argument('resdir', type=str, help="Directory of results")
    psr.add_argument('libpath', type=str, help="Path to library")
    psr.add_argument('num', type=int, help="Number of spectra used in lincomb")
    psr.add_argument('-s', '--suffix', type=str, default="", help="Suffix on results file")
    args = psr.parse_args()

    respath = os.path.join(args.resdir+args.name+"/"+args.name+args.suffix+"_lincomb.csv")
    outpath = os.path.join(args.resdir+args.name+"/"+args.name+args.suffix+"_lincomb{0:d}.pdf".format(args.num))

    main(args.name, respath, outpath, args.libpath, args.num)
