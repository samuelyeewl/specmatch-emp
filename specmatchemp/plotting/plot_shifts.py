"""
@filename plot_shifts.py

Plot shift results
"""

import h5py
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


import os
from argparse import ArgumentParser

from specmatchemp.plotting import plots
from specmatchemp.io import specmatchio

REGION1 = (5158,5172)
REGION2 = (5846,5860)

def main(specpath, outpath, nso=None):
    f = h5py.File(specpath, 'r')
    
    with PdfPages(outpath) as pdf:
        # Plot different regions of the spectrum
        fig = plt.figure(figsize=(10,7))
        gs = gridspec.GridSpec(2,1)
        plt.subplot(gs[0])
        plot_shifts(f, REGION1, nso)
        plt.subplot(gs[1])
        plot_shifts(f, REGION2, nso)
        plt.suptitle('Shift results for star {0}'.format(f.attrs['cps_name']), fontsize=16)
        gs.tight_layout(fig, rect=[0, 0, 1, 0.95])
        pdf.savefig()
        plt.close()

        # Plot lags
        fig = plt.figure(figsize=(8,6))
        plt.title('Lags for shifting of star {0}, spectrum {1}'.format(f.attrs['cps_name'], f.attrs['obs']), fontsize=16)
        plot_lags(f)
        plt.tight_layout()
        pdf.savefig()
        plt.close()

def plot_shifts(f, wavlim, nso):
    """Wrapper function to take in a file object
    """
    s = f['s'][:]
    w = f['w'][:]
    w, s = specmatchio.truncate_spectrum(wavlim, w, s)
    s_un = np.reshape(f['s_unshifted'][:],-1)
    w_un = np.reshape(f['w_unshifted'][:],-1)
    w_un, s_un=  specmatchio.truncate_spectrum(wavlim, w_un, s_un)
    s_ref = f['s_ref'][:]
    w_ref = f['w'][:]
    w_ref, s_ref = specmatchio.truncate_spectrum(wavlim, w_ref, s_ref)
    if nso is not None:
        w_nso, s_nso, serr_nso, hdr_nso = specmatchio.read_standard_spectrum(nso,wavlim=wavlim)
    else:
        s_nso = None
        w_nso = None

    plots.plot_shifts(s, w, s_un, w_un, s_ref, w_ref, s_nso, w_nso, \
        labels={'targ_label':f.attrs['obs'][:], 'ref_label':f.attrs['ref'], 'nso_label':''})

def plot_lags(f):
    plots.plot_lags(f['lag'][:], f['center_pix'][:], f['fit'][:])



if __name__ == '__main__':
    psr = ArgumentParser(description="Build the SpecMatch-Emp library from the various catalogs")
    psr.add_argument('name', type=str, help="CPS name of star")
    psr.add_argument('resdir', type=str, help="Directory of results")
    psr.add_argument('-n', '--nso', type=str, default=None, help="Path to NSO spectrum")
    psr.add_argument('-s', '--suffix', type=str, default="", help="Suffix on results file")
    args = psr.parse_args()

    specpath = os.path.join(args.resdir+args.name+"/"+args.name+args.suffix+"_spec.h5")
    outpath = os.path.join(args.resdir+args.name+"/"+args.name+args.suffix+"_shifts.pdf")
    main(specpath, outpath, args.nso)

