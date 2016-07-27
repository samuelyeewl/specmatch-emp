#!/usr/bin/env python
"""
@filename generate_shift_script.py

Generate script lines for shifting spectra
"""

import os
import re
import pandas as pd

EXECPATH = '/home/syee/specmatchemp-working/specmatchemp/buildlib/shift_spectrum.py'
LIBPATH = '/home/syee/specmatchemp-working/specmatchemp/lib/libstars.csv'
SPECDIR = '/home/syee/specmatchemp-working/specmatchemp/spectra/iodfitsdb/'
REFDIR = '/home/syee/specmatchemp-working/specmatchemp/spectra/refs/'
SHIFTINSTRUCTIONS = '/home/syee/specmatchemp-working/specmatchemp/spectra/shift_reference.csv'
OUTDIR = '/home/syee/specmatchemp-working/specmatchemp/results/'
SCRIPTPATH = '/home/syee/specmatchemp-working/specmatchemp/buildlib/shift_script.txt'


if __name__ == '__main__':
    libstars = pd.read_csv(LIBPATH, index_col=0)
    shiftinst = pd.read_csv(SHIFTINSTRUCTIONS)

    f = open(SCRIPTPATH, 'w')

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
            s+= OUTDIR
            s+= "\n"

            f.write(s)

    libstars.to_csv(LIBPATH)

    f.close()


