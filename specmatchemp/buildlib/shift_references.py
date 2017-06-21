#!/usr/bin/env python
"""
@filename shift_references.py

Shifts spectra onto a common log-lambda scale
"""

import os
from shutil import copyfile

from specmatchemp import SPECMATCHDIR
from specmatchemp import spectrum
from specmatchemp import shift

# List of reference spectra, cps name, teff, and previous reference
from specmatchemp import SHIFT_REFERENCES as REFERENCES


def main():
    # Shift the references onto the NSO wavlength scale
    spectra_dir = os.path.join(SPECMATCHDIR, 'spectra/')
    shifted_dir = os.path.join(SPECMATCHDIR, 'shifted_spectra/')

    if not os.path.exists(shifted_dir):
        os.mkdir(shifted_dir)

    # Save base
    r = REFERENCES[0]
    infile = os.path.join(spectra_dir, r[0]+'.fits')
    base_ref = spectrum.read_fits(infile)
    base_w = base_ref.w
    outfile = os.path.join(shifted_dir, r[0]+'_adj.fits')
    copyfile(infile, outfile)

    # Shift each spectrum in order
    for r in REFERENCES[1:]:
        print("Shifting spectrum {0} onto reference {1}".format(
            r[0], r[3]))

        reffile = os.path.join(shifted_dir, r[3]+'_adj.fits')
        ref_spec = spectrum.read_fits(reffile)

        maskfile = os.path.join(SPECMATCHDIR, 'hires_telluric_mask.csv')

        # Shift each chip
        chips = ['b', 'r', 'i']
        shifted = []
        for c in chips:
            infile = os.path.join(spectra_dir, c + r[0] + '.fits')
            outfile = os.path.join(shifted_dir, c + r[0] + '_adj.fits')
            targ_unshifted = spectrum.read_hires_fits(infile, maskfile=maskfile)

            shift_data = {}

            targ_shifted = shift.shift(targ_unshifted, ref_spec, store=shift_data)

            shift.save_shift_to_fits(outfile, targ_shifted, targ_unshifted,
                                     shift_data, clobber=True)

            shifted.append(targ_shifted)

        # Flatten spectrum and save
        print("Flattening spectra into single file")
        outfile = os.path.join(shifted_dir, r[0] + '_adj.fits')
        flattened = spectrum.Spectrum.combine_spectra(shifted, base_w)
        flattened.to_fits(outfile)


if __name__ == '__main__':
    main()
