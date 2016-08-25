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
                     default='./shift_script.txt',
                     help="Path to output shift script")
    args = psr.parse_args()

    params = pd.read_csv(args.libpath)

    with open(args.outpath, 'w') as f:
        


EXECPATH = '/home/syee/specmatchemp-working/specmatchemp/buildlib/shift_spectrum.py'
LIBPATH = '/home/syee/specmatchemp-working/specmatchemp/lib/libstars.csv'
SPECDIR = '/home/syee/specmatchemp-working/specmatchemp/spectra/iodfitsdb/'
REFDIR = '/home/syee/specmatchemp-working/specmatchemp/spectra/refs/'
TELLURICMASKPATH = '/home/syee/specmatchemp-working/specmatchemp/lib/telluric_mask.csv'
# SHIFTINSTRUCTIONS = '/home/syee/specmatchemp-working/specmatchemp/spectra/shift_reference.csv'
OUTDIR = '/home/syee/specmatchemp-working/specmatchemp/results/'
# SCRIPTPATH = '/home/syee/specmatchemp-working/specmatchemp/buildlib/shift_script.txt'


if __name__ == '__main__':
    psr = ArgumentParser(description="Generate script for shifting")
    psr.add_argument('shiftinstructions', type=str, help="Path to csv file containing shift instructions")
    psr.add_argument('scriptpath', type=str, help="Path to save script")
    psr.add_argument('-s', '--suffix', type=str, default="", help="Suffix to append to result filename")
    args = psr.parse_args()

    libstars = pd.read_csv(LIBPATH, index_col=0)
    shiftinst = pd.read_csv(args.shiftinstructions)
    
    f = open(args.scriptpath, 'w')

    for i in range(len(shiftinst)):
        refpath = os.path.join(REFDIR, shiftinst.iloc[i].ref_spectrum)
        min_t = shiftinst.iloc[i].min_t
        max_t = shiftinst.iloc[i].max_t
        query = '{0:.0f} <= Teff < {1:.0f}'.format(min_t, max_t)
        lib_cut = libstars.query(query)

        for idx, row in lib_cut.iterrows():
            # check for obs
            obs_list = re.sub("[\[\]\'']", '', row.obs).split()
            specpath = None
            for obs in obs_list:
                specpath = os.path.join(SPECDIR, obs+'.fits')
                if os.path.isfile(specpath):
                    libstars.loc[idx, 'lib_obs'] = obs
                else:
                    specpath = None
            if specpath is None:
                print("Could not find any spectra for star {0}".format(targ_params.cps_name))
                continue

            # write script line
            s = "source ~/.bash_profile; "
            s+= "python "+EXECPATH+" "
            s+= row['cps_name']+" "
            s+= specpath+" "
            s+= refpath+" "
            s+= OUTDIR+" "
            s+= "-m "+TELLURICMASKPATH+" "
            if len(args.suffix) > 0:
                s+= " -s "+args.suffix
            s+= "\n"

            f.write(s)

    libstars.to_csv(LIBPATH)

    f.close()


