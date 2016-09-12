#!/usr/bin/env python
"""
@filename generate_shift_script.py

Generate script lines for shifting spectra
"""

from __future__ import print_function

import os
import pandas as pd
from argparse import ArgumentParser

from specmatchemp import SPECMATCHDIR

if __name__ == '__main__':
    psr = ArgumentParser(description="Generate script for shifting")
    psr.add_argument('-l', '--libpath', type=str,
                     default=os.path.join(SPECMATCHDIR, 'libstars.csv'),
                     help="Path to parameters csv file")
    psr.add_argument('-o', '--outpath', type=str,
                     default='./shift_script.sh',
                     help="Path to output shift script")
    psr.add_argument('-s', '--suffix', type=str, default="",
                     help="Suffix to append to shift results")
    args = psr.parse_args()

    params = pd.read_csv(args.libpath)

    with open(args.outpath, 'w') as f:
        for idx, row in params.iterrows():
            obs = row['lib_obs'][1:]
            name = row['cps_name']
            s = "smemp shift " + obs + " "
            s += "-d ~/Dropbox/SpecMatch-Emp/spectra/iodfitsdb/ "
            s += "-p -o ~/SpecMatch-Emp/results/ "
            s += "-n " + name
            s += "\n"
            f.write(s)
