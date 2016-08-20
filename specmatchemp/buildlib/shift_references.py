#!/usr/bin/env python
"""
@filename shift_references.py

Shifts spectra onto a common log-lambda scale
"""

import os
from shutil import copyfile
from astropy.io import fits

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

    # Shift each spectrum in order
    for r in REFERENCES:
        infile = os.path.join(spectra_dir, r[0]+'.fits')
        outfile = os.path.join(shifted_dir, r[0]+'_adj.fits')

        if r[3] is None:
            # No reference to shift against, use this as base.
            copyfile(infile, outfile)
            continue

        print("Shifting spectrum {0} onto reference {1}".format(
            r[0], r[3]))

        reffile = os.path.join(shifted_dir, r[3]+'_adj.fits')
        ref_spec = spectrum.read_fits(reffile)

        targ_unshifted = spectrum.read_hires_fits(infile)

        shift_data = {}

        targ_shifted = shift.shift(targ_unshifted, ref_spec, store=shift_data)

        # store into fits
        prihdu = fits.PrimaryHDU(header=targ_shifted.header)
        targ_shifted_hdu = targ_shifted.to_hdu()

        # store unshifted spectrum
        targ_unshifted_s = fits.ImageHDU(data=targ_unshifted.s,
                                         header=targ_unshifted.header)
        targ_unshifted_serr = fits.ImageHDU(data=targ_unshifted.serr)
        targ_unshifted_w = fits.ImageHDU(data=targ_unshifted.w)

        if len(shift_data) > 0:
            col_list = []
            num_orders = shift_data.pop('num_orders')
            col_list.append(fits.Column(name='num_orders', format='J',
                                        array=[num_orders]))

            num_sects = []
            for i in range(num_orders):
                num_sects.append(shift_data.pop('order_{0:d}/num_sections'
                                                .format(i)))
            col_list.append(fits.Column(name='num_sects', format='J',
                                        array=num_sects))
            n_sec = max(num_sects)

            for k in ['center_pix', 'lag', 'fit']:
                col_list.append(fits.Column(name=k,
                                format='{0:d}E'.format(n_sec),
                                array=shift_data.pop(k)))

            # Save individual fit data
            for k in shift_data.keys():
                col_list.append(fits.Column(name=k, format='D',
                                            array=shift_data[k]))

            shift_hdu = fits.BinTableHDU.from_columns(col_list)
            shift_hdu.name = 'SHIFTDATA'

        hdulist = fits.HDUList([prihdu, targ_shifted_hdu, shift_hdu,
                                targ_unshifted_s, targ_unshifted_serr,
                                targ_unshifted_w])

        hdulist.writeto(outfile)


if __name__ == '__main__':
    main()
