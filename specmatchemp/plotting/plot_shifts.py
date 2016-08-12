"""
@filename plot_shifts.py

Plot shift results
"""

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


import os
from argparse import ArgumentParser

from specmatchemp import spectrum
from specmatchemp import plots
from specmatchemp.io import specmatchio

REGION1 = (5158,5172)
REGION2 = (5846,5860)
ORDER = 11

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

        # Plot xcorr
        fig = plt.figure(figsize=(8,6))
        plt.title('Cross-correlation array for star {0}, order {1:d}'.format(f.attrs['cps_name'], ORDER), fontsize=16)
        plot_xcorr(f, ORDER)
        plt.tight_layout()
        pdf.savefig()
        plt.close()


def plot_shifts(f, wavlim, nsopath=None):
    """Wrapper function to take in a file object
    """
    target = spectrum.read_hdf(f).cut(*wavlim)
    target.plot(offset=1, text='Target (shifted): {0}'.format(target.attrs['obs']))

    unshifted = spectrum.HiresSpectrum(f['w_unshifted'][:], f['s_unshifted'][:]).cut(*wavlim)
    unshifted.plot(offset=0, normalize=True, text='Target (unshifted)', plt_kw={'color':'forestgreen'})

    reference = spectrum.Spectrum(f['w'][:], f['s_ref'][:], attrs={'obs':f.attrs['ref']}).cut(*wavlim)
    reference.plot(offset=1.5, text='Reference: {0}'.format(reference.attrs['obs']), plt_kw={'color':'firebrick'})

    if nsopath is not None:
        nso = spectrum.read_fits(nsopath).cut(*wavlim)
        nso.plot(offset=2.0, text='NSO', plt_kw={'color':'c'})

    plt.plot(target.w, reference.s-target.s, '-', color='purple')
    plots.annotate_spectrum('Residuals', spec_offset=-1)

    plt.ylim(-0.4, 3.4)


def plot_lags(f):
    lags = f['lag'][:]
    center_pix = f['center_pix'][:]
    fit = f['fit'][:]

    num_orders = lags.shape[0]
    # set different colors for each set
    colormap = plt.cm.nipy_spectral

    for i in range(num_orders):
        plt.plot(center_pix[i], lags[i], 'o', color=colormap(0.9*i/num_orders))
        plt.plot(center_pix[i], fit[i], '-', color=colormap(0.9*i/num_orders), label='{0:d}'.format(i))

    plt.xlabel('Pixel number')
    plt.ylabel('Shift (pixels)')
    plt.legend(loc='best', ncol=2, fontsize='small')

def plot_xcorr(f, order):
    """Plot the correlation array produced when shifting

    Args:
        f: h5py file containing xcorr data
        order: Order on chip to plot
    """
    num_sects = len(f['xcorr/order_{0:d}'.format(order)])
    for i in range(num_sects):
        grp = f['xcorr/order_{0:d}/sect_{1:d}'.format(order, i)]
        plt.plot(grp['lag_arr'][:], grp['xcorr'][:], label="{0:d}".format(i))
        max_corr = np.argmax(grp['xcorr'][:])
        plt.plot(grp['lag_arr'][:][max_corr], grp['xcorr'][:][max_corr], 'ko')
    plt.legend(loc='upper left', fontsize='small')



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

