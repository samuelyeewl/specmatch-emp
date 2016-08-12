#!/usr/bin/env python
"""
@filename shift_spectrum.py

Shifts spectra onto a common log-lambda scale
"""

from __future__ import print_function

import os, sys
from argparse import ArgumentParser

from specmatchemp import spectrum
from specmatchemp.shift import shift

import h5py
import numpy as np
import pandas as pd


def main(name, specpath, refpath, outdir, maskpath, suffix):
    targ = spectrum.read_hires_fits(specpath, maskfile=maskpath)
    ref = spectrum.read_fits(refpath)

    # create diagnostic file
    outdir = os.path.join(outdir, name)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    filepath = os.path.join(outdir, name+suffix+'_spec.h5')
    f = h5py.File(filepath, 'w')

    shifted = shift(targ, ref, store=f)
    shifted.to_hdf(f)
    # get wavelength limits
    w_min = shifted.w[0]
    w_max = shifted.w[-1]
    ref_trunc = ref.cut(w_min, w_max)
    f.create_dataset('s_ref', data=ref_trunc.s)
    f.create_dataset('serr_ref', data=ref_trunc.serr)
    # store unshifted spectrum
    f.create_dataset('s_unshifted', data=targ.s)
    f.create_dataset('serr_unshifted', data=targ.serr)
    f.create_dataset('w_unshifted', data=targ.w)

    # store metadata
    f.attrs['cps_name'] = name
    f.attrs['obs'] = os.path.basename(specpath)
    f.attrs['ref'] = os.path.basename(refpath)

    f.close()

if __name__ == '__main__':
    psr = ArgumentParser(description="Build the SpecMatch-Emp library from the various catalogs")
    psr.add_argument('name', type=str, help="CPS name of star")
    psr.add_argument('specpath', type=str, help="Path to target spectrum")
    psr.add_argument('refpath', type=str, help="Path to reference spectrum")
    psr.add_argument('outdir', type=str, help="Directory to output result")
    psr.add_argument('-m', '--mask', type=str, default="", help="Path to telluric mask")
    psr.add_argument('-s', '--suffix', type=str, default="", help="Suffix to append to result filename")
    args = psr.parse_args()

    if not os.path.isfile(args.specpath):
        print("Could not find "+args.specpath)
        sys.exit(1)
    if not os.path.isfile(args.refpath):
        print("Could not find "+args.refpath)
        sys.exit(1)

    print("Shifting star {0}, obs {1} ref {2}".format(args.name, args.specpath, args.refpath))

    main(args.name, args.specpath, args.refpath, args.outdir, args.mask, args.suffix)