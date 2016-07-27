#!/usr/bin/env python
"""
@filename combine_library.py

Combines the library parameters and spectra into a library.h5 file
"""

import os, sys
from argparse import ArgumentParser

import h5py
import numpy as np
import pandas as pd

from specmatchemp import library
from specmatchemp.io import specmatchio

WAVLIM = (4980, 6420)

def main(parampath, specdir, outpath, maskpath):
    libparams = pd.read_csv(parampath, index_col=0)
    wav = None
    spectra = None

    for idx, row in libparams.iterrows():
        # get shifted spectrum
        specpath = os.path.join(specdir,row.cps_name+"/"+row.cps_name+".h5")
        f = h5py.File(specpath, 'r')
        s = f['s'][:]
        serr = f['serr'][:]
        w = f['w'][:]
        f.close()
        # truncate spectrum
        w, s, serr = specmatchio.truncate_spectrum(WAVLIM, w, s, serr)
        if wav is None:
            wav = w
            spectra = np.empty(0,2,len(wav))
        else:
            assert np.allclose(wav, w), "Library spectra not on same wavelength scale"

        # add spectrum to library
        libparams.loc[idx, 'lib_index'] = len(spectra)
        spectra = np.vstack((spectra, [[s, serr]]))

        # calculate signal to noise
        libparams.loc[idx, 'snr'] = np.nanpercentile(1/serr, 90)

    mask = None
    if maskpath is not None:
        mask = pd.read_csv(maskpath, index_col=0)

    # save as library object
    lib = library.Library(wav, spectra, libparams, wavlim=WAVLIM, param_mask=mask)
    lib.to_hdf(outpath)



if __name__ == '__main__':
    psr = ArgumentParser(description="Combined the library parameters and spectra into a library object")
    psr.add_argument('parampath', type=str, help="Path to parameter file")
    psr.add_argument('specdir', type=str, help="Directory containing shifted spectra")
    psr.add_argument('outpath', type=str, help="Path to save library")
    psr.add_argument('-m', '--maskpath', type=str, default=None, help="Path to mask file")
    args = psr.parse_args()

    main(args.parampath, args.specdir, args.outpath, args.maskpath)