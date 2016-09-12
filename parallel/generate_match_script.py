#!/usr/bin/env python
"""
@filename generate_match_script.py

Generate script lines for shifting spectra
"""

from __future__ import print_function

import os
import pandas as pd
from argparse import ArgumentParser

from specmatchemp import library
from specmatchemp import SPECMATCHDIR

if __name__ == '__main__':
    psr = ArgumentParser(description="Generate script for library match")
    psr.add_argument('-l', '--libpath', type=str,
                     default=os.path.join(SPECMATCHDIR, 'library.h5'),
                     help="Path to parameters csv file")
    psr.add_argument('-o', '--outpath', type=str,
                     default='./match_script.sh',
                     help="Path to output match script")
    psr.add_argument('-s', '--suffix', type=str, default="",
                     help="Suffix to append to match results")
    args = psr.parse_args()

    lib = library.read_hdf(args.libpath, wavlim='none')
    params = lib.library_params

    with open(args.outpath, 'w') as f:
        for idx, row in params.iterrows():
            obs = row['lib_obs'][1:]
            name = row['cps_name']
            s = "source ~/.bash_profile; "
            s += "smemp match " + obs + " "
            s += "-pp "
            s += "-i " + name + " "
            s += "-o /home/syee/specmatchemp-working/specmatchemp/results"
            s += "\n"
            f.write(s)
