#!/usr/bin/env python
"""
@filename combine_library.py

Combines the library parameters and spectra into a library.h5 file
"""

import os
from argparse import ArgumentParser

import h5py
import numpy as np
import pandas as pd

from specmatchemp import library
from specmatchemp import spectrum

WAVLIM = (4990, 6410)

def main(parampath, specdir, outpath, maskpath):
    libparams = pd.read_csv(parampath, index_col=0)
    wav = None
    spectra = None

    for idx, row in libparams.iterrows():
        # get shifted spectrum
        specpath = os.path.join(specdir,row.cps_name+"/"+row.cps_name+"_spec.h5")
        spec = spectrum.read_hdf(specpath)
        spec = spec.cut(*WAVLIM)
        if wav is None:
            wav = spec.w
            spectra = np.empty((0,3,len(wav)))
        else:
            assert np.allclose(wav, w), "Library spectra not on same wavelength scale"

        # add spectrum to library
        libparams.loc[idx, 'lib_index'] = len(spectra)
        spectra = np.vstack((spectra, [[spec.s, spec.serr, spec.mask]]))

        # calculate signal to noise
        libparams.loc[idx, 'snr'] = np.nanpercentile(1/spec.serr, 90)

    param_mask = None
    if maskpath is not None:
        param_mask = pd.read_csv(maskpath, index_col=0)

    libparams.drop('obs', axis=1, inplace=True)

    # save as library object
    lib = library.Library(wav, spectra, libparams, wavlim=WAVLIM, param_mask=param_mask)
    lib.to_hdf(outpath)



if __name__ == '__main__':
    psr = ArgumentParser(description="Combined the library parameters and spectra into a library object")
    psr.add_argument('parampath', type=str, help="Path to parameter file")
    psr.add_argument('specdir', type=str, help="Directory containing shifted spectra")
    psr.add_argument('outpath', type=str, help="Path to save library")
    psr.add_argument('-m', '--maskpath', type=str, default=None, help="Path to mask file")
    args = psr.parse_args()

    main(args.parampath, args.specdir, args.outpath, args.maskpath)