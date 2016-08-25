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

from specmatchemp import SPECMATCHDIR
from specmatchemp import SHIFT_REFERENCES
from specmatchemp import spectrum
from specmatchemp import shift
from specmatchemp import plots
from specmatchemp import specmatch
from specmatchemp import library


def match_spectrum(args):
    if not os.path.exists(args.spectrum):
        raise ValueError(args.spectrum + " does not exist!")

    target = spectrum.read_hires_fits(args.spectrum)
    lib = library.read_hdf()
    sm = specmatch.SpecMatch(target, lib)
    sm.shift()
    sm.match()
    sm.lincomb()

    # Print results
    print("SpecMatch Results for {0}".format(os.path.basename(args.spectrum)))
    for p in library.Library.STAR_PROPS:
        print("{0}: {1:.2f}".format(p, sm.results[p]))

    if len(args.name) > 0:
        name = args.name
    else:
        name = os.path.basename(args.spectrum)[:-4]

    # Save results
    outpath = name+'_sm.hdf'
    sm.to_hdf(outpath)

    # Plot results
    if args.plots:
        plotspath = name+'_plots.pdf'
        with PdfPages('./'+plotspath) as pdf:
            # ----------- Shift Plots ---------- #
            wavlim = (5160, 5190)
            order = 2
            # Shifted spectrum
            fig = plt.figure(figsize=(10, 5))
            sm.plot_shifted_spectrum(wavlim=wavlim)
            plt.title('{0} shift results'.format(name))
            fig.set_tight_layout(True)
            pdf.savefig()
            plt.close()

            fig = plt.figure(figsize=(10, 5))
            sm.plot_xcorr(order, True)
            plt.title('{0} cross-correlation for order {0:d}'
                      .format(name, order))
            fig.set_tight_layout(True)
            pdf.savefig()
            plt.close()

            fig = plt.figure(figsize=(8, 6))
            sm.plot_shift_lags()
            plt.title('{0} lags'.format(name))
            fig.set_tight_layout(True)
            pdf.savefig()
            plt.close()

            # ------------ Match plots ------------
            fig = plt.figure(figsize=(12, 8))
            sm.plot_chi_squared_surface()
            plt.title('{0} chi-squared surface'.format(name))
            fig.set_tight_layout(True)
            pdf.savefig()
            plt.close()

            fig = plt.figure(figsize=(10, 4))
            sm.plot_best_match_spectra(region=(5100, 5200),
                                       wavlim=(5160, 5190), num_best=1)
            plt.title('{0} best matching spectrum'.format(name))
            plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])
            pdf.savefig()
            plt.close()

            # ------------ Lincomb plots ------------
            fig = plt.figure(figsize=(10, 8))
            sm.plot_references(region=(5100, 5200), verbose=True)
            axes = fig.axes
            axes[0].legend(numpoints=1, fontsize='small', loc='best')
            plt.title('{0} references used in linear combination'.format(name))
            plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])
            pdf.savefig()
            plt.close()

            fig = plt.figure(figsize=(12, 6))
            sm.plot_lincomb(region=(5100, 5200), wavlim=(5160, 5190))
            plt.title('{0} Linear Combination results'.format(name))
            fig.set_tight_layout(True)
            pdf.savefig()
            plt.close()

    return sm


def shift_spectrum(args):
    """Shift a target spectrum given an observation code.

    Saves the shifted spectrum in a fits file.

    Returns:
        shifted, unshifted, shift_data
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
    outpath = os.path.join(shiftedspecdir, 'r'+args.obs+'_adj.fits')
    shift.save_shift_to_fits(outpath, shifted, targ_spec, shift_data)

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

    return shifted, targ_spec, shift_data


def main():
    args = sys.argv[1:]

    psr = ArgumentParser(description="SpecMatch-Emp: A tool for extracting " +
                         "fundamental stellar parameters from spectra, updated",
                         prog='smemp')
    subpsr = psr.add_subparsers(title="subcommands", dest='subcommand',
                description="specmatch: Perform the entire SpecMatch " +
                            "algorithm on a HIRES spectrum" +
                            "shift: Shift a spectrum onto the reference " +
                            "wavelength scale.")
    subpsr.required = True

    psr_sm = subpsr.add_parser("specmatch")
    psr_sm.add_argument("spectrum", type=str, help="Path to spectrum file")
    psr_sm.add_argument("-p", "--plots", action='store_true',
                        help="Generate diagnostic plots")
    psr_sm.add_argument("-n", "--name", type=str, default="",
                        help="Star name")
    psr_sm.set_defaults(func=match_spectrum)

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
