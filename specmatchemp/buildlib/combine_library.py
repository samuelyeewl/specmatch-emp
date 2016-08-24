#!/usr/bin/env python
"""
@filename combine_library.py

Combines the library parameters and spectra into a library.h5 file
"""

from __future__ import print_function

import os
from argparse import ArgumentParser

import h5py
import numpy as np
import pandas as pd

from specmatchemp import SPECMATCHDIR
from specmatchemp import SHIFT_REFERENCES
from specmatchemp import library
from specmatchemp import spectrum

WAVLIM = (4990, 6410)


def main(paramdir, specdir, outdir):
    library_params = pd.read_csv(os.path.join(paramdir, 'libstars.csv'),
                                 index_col=0)

    library_params['lib_index'] = None
    library_params['lib_obs'] = None

    # Read in NSO spectrum
    nso = spectrum.read_fits(os.path.join(specdir, 'nso_adj.fits'))
    nso = nso.cut(*WAVLIM)
    # Use NSO wavelenght scale as library wavelength scale
    wav = nso.w

    # Read in all spectra
    spectra = np.empty((0, 3, len(wav)))
    for idx, row in library_params.iterrows():
        lib_obs = None
        # Check for all observations
        for obs in row.obs:
            if os.path.exists(os.path.join(specdir, obs+'_adj.fits')):
                lib_obs = obs
        if lib_obs is None:
            print("Could not find spectrum for star {0}".format(row.cps_name))
            continue

        library_params.loc[idx, 'lib_obs'] = lib_obs
        spec_path = os.path.join(specdir, lib_obs+'_adj.fits')
        spec = spectrum.read_fits(spec_path).cut(*WAVLIM)

        # Check that shifted spectra are on same wavelength scale
        if not np.allclose(spec.w, wav):
            print("Spectrum {0} was not on the same wavelength scale".format(
                  row.lib_obs))
            library_params.loc[idx, 'lib_obs'] = None
            continue

        # Add spectrum to library
        library_params.loc[idx, 'lib_index'] = len(spectra)
        spectra = np.vstack((spectra, [[spec.s, spec.serr, spec.mask]]))

        # Calculate signal to noise
        library_params.loc[idx, 'snr'] = np.nanpercentile(1/spec.serr, 90)

    # Read in parameter mask
    param_mask = pd.read_csv(os.path.join(paramdir, 'libstars_mask.csv'))

    # Drop stars with missing spectra
    missing_spectra = pd.isnull(library_params.lib_obs)
    missing_indices = library_params[missing_spectra].index.tolist()
    library_params.drop(missing_indices, inplace=True)
    param_mask.drop(missing_indices, inplace=True)

    lib = library.Library(wav, spectra, library_params, wavlim=WAVLIM,
                          param_mask=param_mask, nso=nso)

    # Get allowed shift references
    shift_obs = [row[0] for r in SHIFT_REFERENCES]
    lib.header['shift_refs'] = shift_obs

    outpath = os.path.join(outdir, 'library.h5')
    lib.to_hdf(outpath)


if __name__ == '__main__':
    psr = ArgumentParser(description="Combine the library parameters and " +
                         "spectra into a library object.")
    psr.add_argument('-p', '--paramdir', type=str, default=SPECMATCHDIR,
                     help="Directory to check for parameter file")
    psr.add_argument('-s', '--specdir', type=str,
                     default=os.path.join(SPECMATCHDIR, 'shifted_spectra/'),
                     help="Directory to check for shifted spectra")
    psr.add_argument('-o', '--outdir', type=str, default=SPECMATCHDIR,
                     help="Directory to save output")
    args = psr.parse_args()

    main(args.paramdir, args.specdir, args.outdir)
