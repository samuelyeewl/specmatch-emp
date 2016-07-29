"""
@filename plot_shifts.py

Plot shift results
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import os
from argparse import ArgumentParser

from specmatchemp.plotting import plots
from specmatchemp.io import specmatchio

def main(specpath, min_w, max_w, nso=None):
    wavlim = (min_w, max_w)

    f = h5py.File(specpath, 'r')
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

    plt.figure(figsize=(12,6))
    plots.plot_shifts(s, w, s_un, w_un, s_ref, w_ref, s_nso, w_nso)
    
    # hide y ticks
    plt.title('Shift results for star {0}'.format(f.attrs['cps_name']))
    plt.show()



if __name__ == '__main__':
    psr = ArgumentParser(description="Build the SpecMatch-Emp library from the various catalogs")
    psr.add_argument('name', type=str, help="CPS name of star")
    psr.add_argument('resdir', type=str, help="Directory of results")
    psr.add_argument('min_w', type=float, help="Minimum wavelength to plot")
    psr.add_argument('max_w', type=float, help="Maximum wavelength to plot")
    psr.add_argument('-n', '--nso', type=str, default=None, help="Path to NSO spectrum")
    psr.add_argument('-s', '--suffix', type=str, default="", help="Suffix on results file")
    args = psr.parse_args()

    specpath = os.path.join(args.resdir+args.name+"/"+args.name+args.suffix+"_spec.h5")
    main(specpath, args.min_w, args.max_w, args.nso)

