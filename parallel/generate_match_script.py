#!/usr/bin/env python

from argparse import ArgumentParser
from specmatchemp import library

WAV_MIN = 5000
WAV_MAX = 6400
WAV_STEP = 100

LIBPATH = '/home/syee/specmatchemp-working/specmatchemp/lib/library.h5'
EXECPATH = '/home/syee/specmatchemp-working/specmatchemp/tests/match_library_spectrum.py'
OUTDIR = '/home/syee/specmatchemp-working/specmatchemp/results/'

if __name__ == '__main__':
    psr = ArgumentParser(description="Generate script for parallelization")
    psr.add_argument('libpath', type=str, help="Path to library file")
    psr.add_argument('outpath', type=str, help="Path to output file")
    args = psr.parse_args()

    lib = library.read_hdf(args.libpath, wavlim='none')

    f = open(args.outpath, "w")

    for wl in range(WAV_MIN, WAV_MAX, WAV_STEP):
        for name in lib.library_params.cps_name:
            s = "source ~/.bash_profile; "
            s+= "python "+EXECPATH+" "
            s+= LIBPATH+" "
            s+= name+" "
            s+= OUTDIR+" "
            s+= "{0:d} {1:d} ".format(wl, WAV_STEP)
            s+= "\n"

            f.write(s)


    f.close()