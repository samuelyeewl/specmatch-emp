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


def main(name, specpath, refpath, outdir, suffix):
    w_targ, s_targ, serr_targ, hdr_targ = specmatchio.read_hires_spectrum(specpath)
    w_ref, s_ref, serr_ref, hdr_ref = specmatchio.read_standard_spectrum(refpath)

    # create diagnostic file
    outdir = os.path.join(outdir, '/'+name)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    filepath = os.path.join(outdir, name+'.h5')
    f = h5py.File(filepath, 'w')

    
    


if __name__ == '__main__':
    psr = ArgumentParser(description="Build the SpecMatch-Emp library from the various catalogs")
    psr.add_argument('name', type=str, help="CPS name of star")
    psr.add_argument('specpath', type=str, help="Path to target spectrum")
    psr.add_argument('refpath', type=str, help="Path to reference spectrum")
    psr.add_argument('outdir', type=str, help="Directory to output result")
    psr.add_argument('-s', '--suffix', type=str, default="", help="Suffix to append to result filename")
    args = psr.parse_args()

    if not os.isfile(args.specpath):
        print("Could not find "+args.specpath)
        sys.exit(1)
    if not os.isfile(args.refpath):
        print("Could not find "+args.refpath)
        sys.exit(1)

    main(args.name, args.specpath, args.refpath, args.outdir, args.suffix)