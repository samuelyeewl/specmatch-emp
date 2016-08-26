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


def specmatch_spectrum(args):
    if not os.path.exists(args.spectrum):
        raise ValueError(args.spectrum + " does not exist!")

    target = spectrum.read_hires_fits(args.spectrum)
    lib = library.read_hdf()
    sm = specmatch.SpecMatch(target, lib)
    sm.shift()

    if args.in_library:
        sm.target.name = args.in_library
        name = args.in_library
        targ_idx = lib.get_index(args.in_library)
        targ_param, targ_spec = lib[targ_idx]
        sm.match(ignore=targ_idx)
    else:
        name = os.path.basename(args.spectrum)[:-4]
        sm.match()
    sm.lincomb(args.num_best)

    # Print results
    print("SpecMatch Results for {0}".format(name))
    for p in library.Library.STAR_PROPS:
        print("{0}: {1:.2f}".format(p, sm.results[p]))

    outdir = os.path.join(args.outdir, name)
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # Save final results
    outpath = os.path.join(outdir, name + args.suffix + '_results.txt')
    with open(outpath, 'w') as f:
        f.write('Derived Parameters\n')
        f.write('------------------\n')
        f.write('Teff: {0:.0f} +/- {1:.0f} K\n'.format(
            sm.results['Teff'], sm.results['u_Teff']))
        f.write('Radius: {0:.3f} +/- {1:.3f} Rsun\n'.format(
            sm.results['radius'], sm.results['u_radius']))
        f.write('[Fe/H]: {0:.2f} +/- {1:.2f} dex\n'.format(
            sm.results['feh'], sm.results['u_feh']))

        f.write('\n')
        if args.in_library:
            f.write('Library Parameters\n')
            f.write('------------------\n')
            f.write('Teff: {0:.0f} +/- {1:.0f} K\n'.format(
                targ_param['Teff'], targ_param['u_Teff']))
            f.write('Radius: {0:.3f} +/- {1:.3f} Rsun\n'.format(
                targ_param['radius'], targ_param['u_radius']))
            f.write('[Fe/H]: {0:.2f} +/- {1:.2f} dex\n'.format(
                targ_param['feh'], targ_param['u_feh']))

        f.write('\n')
        f.write('Best Matching Spectra\n')
        f.write('---------------------\n')
        for i in range(len(sm.regions)):
            f.write('Region {0}:\n'.format(sm.regions[i]))
            mt = sm.lincomb_matches[i]
            for j in range(mt.num_refs):
                ref = mt.refs[j]
                f.write('\t#{0:d}: {1}, '.format(j, ref.name))
                f.write('chi^2 = {0:.3f}, '.format(mt.ref_chisq[j]))
                f.write('c_{0:d} = {1:.3f}\n'.format(j, mt.coeffs[j]))
            f.write('Final chi^2 = {0:.3f}'.format(mt.best_chisq))

    # Save full results
    outpath = os.path.join(outdir, name + args.suffix + '_sm.hdf')
    sm.to_hdf(outpath)

    # Plot results
    if args.plots:
        plotspath = os.path.join(outdir, name + args.suffix + '_plots.pdf')
        with PdfPages(plotspath) as pdf:
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
            plt.title('{0} cross-correlation for order {1:d}'
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
            if args.in_library:
                # Plot target parameters if available
                axes = fig.axes
                axes[0].axvline(targ_param['Teff'], color='k')
                axes[1].axvline(targ_param['radius'], color='k')
                axes[2].axvline(targ_param['feh'], color='k')
            plt.title('{0} chi-squared surface'.format(name))
            fig.set_tight_layout(True)
            pdf.savefig()
            plt.close()

            fig = plt.figure(figsize=(10, 4))
            sm.plot_best_match_spectra(region=(5100, 5200),
                                       wavlim=wavlim, num_best=1)
            plt.title('{0} best matching spectrum'.format(name))
            plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])
            pdf.savefig()
            plt.close()

            # ------------ Lincomb plots ------------
            fig = plt.figure(figsize=(10, 8))
            sm.plot_references(region=(5100, 5200), verbose=True)
            axes = fig.axes
            axes[0].legend(numpoints=1, fontsize='small', loc='best')
            if args.in_library:
                # Plot target parameters if available
                axes[0].plot(targ_param['Teff'], targ_param['radius'], '*',
                             ms=15, color='red', label='Target')
                axes[1].plot(targ_param['Teff'], targ_param['radius'], '*',
                             ms=15, color='red')
                axes[2].plot(targ_param['feh'], targ_param['radius'], '*',
                             ms=15, color='red')
                axes[3].plot(targ_param['feh'], targ_param['radius'], '*',
                             ms=15, color='red')
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


def match_spectrum(args):
    inpath = os.path.join(args.directory, 'r' + args.obs + '_adj' +
                          args.suffix + '.fits')
    if not os.path.exists(inpath):
        raise ValueError(inpath + " does not exist!")

    # Load shifted spectrum
    target = spectrum.read_fits(inpath)

    lib = library.read_hdf()
    sm = specmatch.SpecMatch(target, lib)

    if args.in_library:
        name = args.in_library
        sm.target.name = args.in_library
        targ_idx = lib.get_index(args.in_library)
        targ_param, targ_spec = lib[targ_idx]
        sm.match(ignore=targ_idx)
    else:
        name = os.path.basename(args.spectrum)[:-4]
        sm.match()

    # Save results
    outdir = os.path.join(args.outdir, name)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    outpath = os.path.join(outdir, name + args.suffix + '_match.csv')
    sm.match_results.to_csv(outpath)

    # Save SpecMatch object
    outpath = os.path.join(outdir, name + args.suffix + '_sm.hdf')
    sm.to_hdf(outpath)

    # Generate plots
    if args.plots:
        plotspath = os.path.join(outdir, name + args.suffix +
                                 '_match_plots.pdf')
        with PdfPages(plotspath) as pdf:
            wavlim = (5160, 5190)
            # ------------ Match plots ------------
            fig = plt.figure(figsize=(12, 8))
            sm.plot_chi_squared_surface()
            if args.in_library:
                # Plot target parameters if available
                axes = fig.axes
                axes[0].axvline(targ_param['Teff'], color='k')
                axes[1].axvline(targ_param['radius'], color='k')
                axes[2].axvline(targ_param['feh'], color='k')
            plt.title('{0} chi-squared surface'.format(name))
            fig.set_tight_layout(True)
            pdf.savefig()
            plt.close()

            fig = plt.figure(figsize=(10, 4))
            sm.plot_best_match_spectra(region=(5100, 5200),
                                       wavlim=wavlim, num_best=1)
            plt.title('{0} best matching spectrum'.format(name))
            plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])
            pdf.savefig()
            plt.close()

    return sm


def lincomb_spectrum(args):
    lib = library.read_hdf()
    sm = specmatch.SpecMatch.read_hdf(args.match_results, lib)

    if args.in_library:
        name = args.in_library
        sm.target.name = args.in_library
        targ_idx = lib.get_index(args.in_library)
        targ_param, targ_spec = lib[targ_idx]
    else:
        name = os.path.basename(args.spectrum)[:-4]

    sm.lincomb(num_best=args.num_best)

    # Print results
    print("SpecMatch Results for {0}".format(name))
    for p in library.Library.STAR_PROPS:
        print("{0}: {1:.2f}".format(p, sm.results[p]))

    outdir = os.path.join(args.outdir, name)
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # Save final results
    outpath = os.path.join(outdir, name + args.suffix + '_results.txt')
    with open(outpath, 'w') as f:
        f.write('Derived Parameters\n')
        f.write('------------------\n')
        f.write('Teff: {0:.0f} +/- {1:.0f} K\n'.format(
            sm.results['Teff'], sm.results['u_Teff']))
        f.write('Radius: {0:.3f} +/- {1:.3f} Rsun\n'.format(
            sm.results['radius'], sm.results['u_radius']))
        f.write('[Fe/H]: {0:.2f} +/- {1:.2f} dex\n'.format(
            sm.results['feh'], sm.results['u_feh']))

        f.write('\n')
        if args.in_library:
            f.write('Library Parameters\n')
            f.write('------------------\n')
            f.write('Teff: {0:.0f} +/- {1:.0f} K\n'.format(
                targ_param['Teff'], targ_param['u_Teff']))
            f.write('Radius: {0:.3f} +/- {1:.3f} Rsun\n'.format(
                targ_param['radius'], targ_param['u_radius']))
            f.write('[Fe/H]: {0:.2f} +/- {1:.2f} dex\n'.format(
                targ_param['feh'], targ_param['u_feh']))

        f.write('\n')
        f.write('Best Matching Spectra\n')
        f.write('---------------------\n')
        for i in range(len(sm.regions)):
            f.write('Region {0:d}:\n'.format(i))
            mt = sm.lincomb_matches[i]
            for j in range(mt.num_refs):
                ref = mt.refs[j]
                f.write('\t#{0:d}: {1}, '.format(j, ref.name))
                f.write(r'$\chi^2 = {0:.3f}$, '.format(mt.ref_chisq[j]))
                f.write(r'$c_{0:d} = {1:.3f}$\n'.format(j, mt.coeffs[j]))
            f.write(r'Final $\chi^2 = {0:.3f}$\n'.format(mt.best_chisq))

    # Save full results
    outpath = os.path.join(outdir, name + args.suffix + '_lincomb_sm.hdf')
    sm.to_hdf(outpath)

    # Plot results
    if args.plots:
        plotspath = os.path.join(outdir, name + args.suffix +
                                 '_lincomb_plots.pdf')
        with PdfPages(plotspath) as pdf:
            # ------------ Lincomb plots ------------
            fig = plt.figure(figsize=(10, 8))
            sm.plot_references(region=(5100, 5200), verbose=True)
            axes = fig.axes
            axes[0].legend(numpoints=1, fontsize='small', loc='best')
            if args.in_library:
                # Plot target parameters if available
                axes[0].plot(targ_param['Teff'], targ_param['radius'], '*',
                             ms=15, color='red', label='Target')
                axes[1].plot(targ_param['Teff'], targ_param['radius'], '*',
                             ms=15, color='red')
                axes[2].plot(targ_param['feh'], targ_param['radius'], '*',
                             ms=15, color='red')
                axes[3].plot(targ_param['feh'], targ_param['radius'], '*',
                             ms=15, color='red')
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
    if args.directory is not specdir:
        copy(os.path.join(args.directory, 'r'+args.obs+'.fits'), specdir)

    # load target and references
    targ_path = os.path.join(specdir, 'r'+args.obs+'.fits')
    targ_spec = spectrum.read_hires_fits(targ_path)

    ref_specs = [spectrum.read_fits(os.path.join(SPECMATCHDIR,
                 'shifted_spectra/'+r[0]+'_adj.fits'))
                 for r in SHIFT_REFERENCES]

    # Shift spectrum
    shift_data = {}
    shifted = shift.bootstrap_shift(targ_spec, ref_specs, store=shift_data)

    # Save shifted spectrum
    outpath = os.path.join(shiftedspecdir, 'r' + args.obs + '_adj' +
                           args.suffix + '.fits')
    shift.save_shift_to_fits(outpath, shifted, targ_spec, shift_data)

    # Generate plots
    if args.plots:
        plotfile = "./" + args.obs + "_shift_plots.pdf"
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
                            "wavelength scale." +
                            "match: Perform only the library grid search " +
                            "portion of the SpecMatch algorithm" +
                            "lincomb: Perform only the linear combination " +
                            "portion of the SpecMatch algorithm")
    subpsr.required = True

    psr_sm = subpsr.add_parser("specmatch")
    psr_sm.add_argument("spectrum", type=str, help="Path to spectrum file")
    psr_sm.add_argument("-p", "--plots", action='store_true',
                        help="Generate diagnostic plots")
    psr_sm.add_argument("-i", "--in_library", type=str, default=False,
                        help="Check the library for a star with the same " +
                        "name and exclude it from the match.")
    psr_sm.add_argument("-o", "--outdir", type=str, default="./",
                        help="Directory to store output files.")
    psr_sm.add_argument("-n", "--num_best", type=int, default=5,
                        help="Number of best matches to use in linear " +
                        "combination step.")
    psr_sm.add_argument("-s", "--suffix", type=str, default="",
                        help="Suffix to append to results files")
    psr_sm.set_defaults(func=specmatch_spectrum)

    psr_shift = subpsr.add_parser("shift")
    psr_shift.add_argument("obs", type=str, help="cps id of target spectrum")
    psr_shift.add_argument("-d", "--directory", type=str,
                           default=os.path.join(SPECMATCHDIR, 'spectra'),
                           help="Directory to look in for spectra")
    psr_shift.add_argument("-p", "--plots", action='store_true',
                           help="Generate diagnostic plots")
    psr_shift.add_argument("-s", "--suffix", type=str, default="",
                        help="Suffix to append to results files")
    psr_shift.set_defaults(func=shift_spectrum)

    psr_match = subpsr.add_parser("match")
    psr_match.add_argument("obs", type=str, help="cps id of target spectrum")
    psr_match.add_argument("-d", "--directory", type=str,
                           default=os.path.join(SPECMATCHDIR, 'shifted_spectra'),
                           help="Directory to look in for spectra")
    psr_match.add_argument("-p", "--plots", action='store_true',
                           help="Generate diagnostic plots")
    psr_match.add_argument("-i", "--in_library", type=str, default=False,
                           help="Check the library for a star with the same " +
                           "name and exclude it from the match.")
    psr_match.add_argument("-o", "--outdir", type=str, default="./",
                        help="Directory to store output files.")
    psr_match.add_argument("-s", "--suffix", type=str, default="",
                        help="Suffix to append to results files")
    psr_match.set_defaults(func=match_spectrum)

    psr_lincomb = subpsr.add_parser("lincomb")
    psr_lincomb.add_argument("match_results", type=str, default="",
                             help="Results hdf file containing match results")
    psr_lincomb.add_argument("-p", "--plots", action='store_true',
                             help="Generate diagnostic plots")
    psr_lincomb.add_argument("-i", "--in_library", type=str, default=False,
                             help="Check the library for a star with the " +
                             "same name and exclude it from the match.")
    psr_lincomb.add_argument("-o", "--outdir", type=str, default="./",
                             help="Directory to store output files.")
    psr_lincomb.add_argument("-s", "--suffix", type=str, default="",
                             help="Suffix to append to results files")
    psr_lincomb.set_defaults(func=lincomb_spectrum)

    args = psr.parse_args(args=args)
    args.func(args)


if __name__ == '__main__':
    main()
