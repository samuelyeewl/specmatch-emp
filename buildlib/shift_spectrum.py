#!/usr/bin/env python
"""
@filename shift_spectrum.py

Shifts spectra onto a common log-lambda scale
"""

import os, sys
from argparse import ArgumentParser

from specmatchemp import shift_spectra
from specmatchemp.io import specmatchio

import h5py
import numpy as np


def main(name, specpath, refpath, outdir, suffix):
    w_targ, s_targ, serr_targ, hdr_targ = specmatchio.read_hires_spectrum(specpath)
    w_ref, s_ref, serr_ref, hdr_ref = specmatchio.read_standard_spectrum(refpath)

    # create diagnostic file
    outdir = os.path.join(outdir, name)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    filepath = os.path.join(outdir, name+'.h5')
    f = h5py.File(filepath, 'w')

    s, serr, w = shift_spectra.shift(s_targ, serr_targ, w_targ, s_ref, serr_ref, w_ref, outfile=f)
    f.create_dataset('s', data=s)
    f.create_dataset('serr', data=serr)
    f.create_dataset('w', data=w)
    # get wavelength limits
    w_min = w[0]
    w_max = w[-1]
    w_ref_trunc, s_ref_trunc, serr_ref_trunc = specmatchio.truncate_spectrum((w_min, w_max), w_ref, s_ref, serr_ref)
    f.create_dataset('s_ref', data=s_ref_trunc)
    f.create_dataset('serr_ref', data=serr_ref_trunc)
    # store unshifted spectrum
    f.create_dataset('s_unshifted', data=s_targ)
    f.create_dataset('serr_unshifted', data=serr_targ)
    f.create_dataset('w_unshifted', data=w_targ)

    f.close()

if __name__ == '__main__':
    psr = ArgumentParser(description="Build the SpecMatch-Emp library from the various catalogs")
    psr.add_argument('name', type=str, help="CPS name of star")
    psr.add_argument('specpath', type=str, help="Path to target spectrum")
    psr.add_argument('refpath', type=str, help="Path to reference spectrum")
    psr.add_argument('outdir', type=str, help="Directory to output result")
    psr.add_argument('-s', '--suffix', type=str, default="", help="Suffix to append to result filename")
    args = psr.parse_args()

    if not os.path.isfile(args.specpath):
        print("Could not find "+args.specpath)
        sys.exit(1)
    if not os.path.isfile(args.refpath):
        print("Could not find "+args.refpath)
        sys.exit(1)

    print("Shifting star {0}, obs {1} ref {2}".format(args.name, args.specpath, args.refpath))

    main(args.name, args.specpath, args.refpath, args.outdir, args.suffix)