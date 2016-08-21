"""
@filename cli.py

Command-line interface
"""
import os
import sys
from shutil import copy
from argparse import ArgumentParser

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

from astropy.io import fits

from specmatchemp import SPECMATCHDIR
from specmatchemp import SHIFT_REFERENCES
from specmatchemp import spectrum
from specmatchemp import shift
from specmatchemp import plots


def shift_spectrum(args):
    """Shift a target spectrum given an observation code.

    Saves the shifted spectrum in a fits file.
    """
    # if a different directory is provided, copy the file into specmatchemp
    # working directory
    specdir = os.path.join(SPECMATCHDIR, 'spectra')
    shiftedspecdir = os.path.join(SPECMATCHDIR, 'shifted_spectra')
    if args.directory is not None:
        copy(os.path.join(args.directory, 'r'+args.obs+'.fits'), specdir)

    # load target and references
    targ_path = os.path.join(specdir, 'r'+args.obs+'.fits')
    targ_spec = spectrum.read_hires_fits(targ_path)

    ref_specs = [spectrum.read_fits(os.path.join(SPECMATCHDIR,
                'shifted_spectra/'+r[0]+'_adj.fits')) for r in SHIFT_REFERENCES]

    # Shift spectrum
    shift_data = {}
    shifted = shift.bootstrap_shift(targ_spec, ref_specs, store=shift_data)

    # Save shifted spectrum
    prihdu = fits.PrimaryHDU(header=shifted.header)
    shifted_hdu = shifted.to_hdu()
    unshifted_hdus = targ_spec.to_hdu()
    shift_data_hdu = shift.shift_data_to_hdu(shift_data)

    hdulist = fits.HDUList([prihdu, shifted_hdu, shift_data_hdu] +
                           unshifted_hdus)

    outpath = os.path.join(shiftedspecdir, 'r'+args.obs+'_adj.fits')
    hdulist.writeto(outpath, clobber=True)

    # Generate plots
    if args.plots:
        plotfile = "./"+args.obs+"_shifts.pdf"
        wavlim = (5158, 5172)

        with PdfPages(plotfile) as pdf:
            targid = targ_spec.name
            fig = plt.figure(figsize=(10, 6))
            plt.title('Shift results for star {0}'.format(targid))
            targ_spec.plot(offset=0, normalize=True, text='Target (unshifted)',
                           plt_kw={'color': 'forestgreen'})

            shifted.plot(offset=1, text='Target (shifted): {0}'.format(targid))

            shift_ref = ref_specs[shift_data['shift_reference']]
            refid = shift_ref.name
            shift_ref.plot(offset=1.5, text='Reference: {0}'.format(refid),
                           plt_kw={'color': 'firebrick'})

            if (shifted.w[0] > shift_ref.w[0]) \
                    or (shifted.w[-1] < shift_ref.w[-1]):
                shifted = shifted.extend(shift_ref.w)
            plt.plot(shifted.w, shift_ref.s - shifted.s, '-', color='purple')
            plots.annotate_spectrum('Residuals', spec_offset=-1)

            plt.xlim(wavlim)
            plt.ylim(-0.5, 2.7)
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close()


def main():
    args = sys.argv[1:]

    psr = ArgumentParser(description="SpecMatch-Emp: A tool for extracting " +
                         "fundamental stellar parameters from spectra",
                         prog='smemp')
    subpsr = psr.add_subparsers(title="subcommands", dest='subcommand',
                description="shift      Shift a spectrum onto the reference " +
                            "wavelength scale.")
    subpsr.required = True

    psr_shift = subpsr.add_parser("shift")
    psr_shift.add_argument("obs", type=str, help="cps id of target spectrum")
    psr_shift.add_argument("-d", "--directory", type=str, default=None,
                           help="Directory to look in for spectra")
    psr_shift.add_argument("-p", "--plots", action='store_true',
                           help="Generate diagnostic plots")
    psr_shift.set_defaults(func=shift_spectrum)

    args = psr.parse_args(args=args)
    args.func(args)


if __name__ == '__main__':
    main()
