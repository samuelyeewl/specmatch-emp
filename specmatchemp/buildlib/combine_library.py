#!/usr/bin/env python
"""
@filename combine_library.py

Combines the library parameters and spectra into a library.h5 file
"""

from __future__ import print_function

import os
from argparse import ArgumentParser

import numpy as np
import pandas as pd

from specmatchemp import SPECMATCHDIR
from specmatchemp import SHIFT_REFERENCES
from specmatchemp import library
from specmatchemp import spectrum

WAVLIM = (3650, 7960)


def main(append, parampath, specdir, outdir):
    library_params = pd.read_csv(parampath, index_col=0)

    library_params['lib_index'] = None

    # Read in NSO spectrum
    nso = spectrum.read_fits(os.path.join(specdir, 'nso_adj.fits'))
    nso = nso.cut(*WAVLIM)
    # Use NSO wavelength scale as library wavelength scale
    wav = nso.w

    # Read in all spectra
    spectra = np.empty((0, 3, len(wav)))
    for idx, row in library_params.iterrows():
        lib_obs = row.lib_obs.lstrip('r')
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

    # Read in parameter mask
    maskpath = parampath[:-4] + '_mask.csv'
    if os.path.exists(maskpath):
        param_mask = pd.read_csv(maskpath, index_col=0)
    else:
        param_mask = pd.DataFrame()
        for p in library.Library.STAR_PROPS:
            param_mask[p] = np.isfinite(library_params[p])
            param_mask['u_'+p] = np.isfinite(library_params[p])

    # Drop stars with missing spectra
    missing_spectra = pd.isnull(library_params.lib_obs)
    missing_indices = library_params[missing_spectra].index.tolist()
    print("Dropping {0:d} stars with no spectra".format(len(missing_indices)))
    library_params.drop(missing_indices, inplace=True)
    param_mask.drop(missing_indices, inplace=True)

    libpath = os.path.join(outdir, 'library.h5')
    if append:
        lib = library.read_hdf(libpath)
        lib.append(library_params, spectra, param_mask)
    else:
        lib = library.Library(wav, spectra, library_params, wavlim=WAVLIM,
                              param_mask=param_mask, nso=nso)
        # Get allowed shift references
        shift_obs = [row[0] for r in SHIFT_REFERENCES]
        lib.header['shift_refs'] = shift_obs

    # Save library
    lib.to_hdf(libpath)


if __name__ == '__main__':
    psr = ArgumentParser(description="Combine the library parameters and " +
                         "spectra into a library object.")
    psr.add_argument('-a', '--append', action='store_true',
                     help="Append to an existing library in outdir.")
    psr.add_argument('-p', '--parampath', type=str,
                     default=os.path.join(SPECMATCHDIR, 'libstars.csv'),
                     help="Directory to check for parameter file")
    psr.add_argument('-s', '--specdir', type=str,
                     default=os.path.join(SPECMATCHDIR, 'shifted_spectra/'),
                     help="Directory to check for shifted spectra")
    psr.add_argument('-o', '--outdir', type=str, default=SPECMATCHDIR,
                     help="Directory to save output")
    args = psr.parse_args()

    main(args.append, args.parampath, args.specdir, args.outdir)
