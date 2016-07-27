#!/usr/bin/env python
"""
@filename generate_shift_script.py

Generate script lines for shifting spectra
"""

import pandas as pd

LIBPATH = '/home/syee/specmatchemp-working/specmatchemp/lib/libstars.csv'
SPECDIR = 

if __name__ == '__main__':
    psr = ArgumentParser(description="Generate script for parallelization")
    psr.add_argument('libpath', type=str, help="Path to library parameters csv file")
    psr.add_argument('specdir', type=str, help="Directory containing spectra")
    psr.add_argument('refdir', type=str, help="Directory containing reference spectra")
    psr.add_argument('shiftinstructions', type=str, help="Path to shift instructions")
    psr.add_argument('outdir', type=str, help="Directory to store output files")
    psr.add_argument('scriptpath', type=str, help="Path to save script")
    args = psr.parse_args()

    libstars = pd.read_csv(args.libpath, index_col=0)
    shiftinst = pd.read_csv(args.shiftinstructions)

    for i in range(len(shiftinst)):
        min_t = shiftinst.iloc[i].min_t
        max_t = shiftinst.iloc[i].max_t
        query = '{0:.0f} <= Teff < {1:.0f}'.format(min_t, max_t)
        lib_cut = libstars.query(query)

        for idx, row in lib_cut.iterrows():
            s = 'source '


