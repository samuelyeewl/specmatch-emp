"""
@filename cli.py

Command-line interface
"""
import os
import sys
from argparse import ArgumentParser

from specmatchemp import core
from specmatchemp import SPECMATCHDIR


def specmatch_spectrum(args):
    core.specmatch_spectrum(args.spectrum, plot_level=args.plots,
                            inlib=args.in_library, outdir=args.outdir,
                            num_best=args.num_best, suffix=args.suffix)


def match_spectrum(args):
    core.match_spectrum(args.spectrum, indir=args.directory,
                        plot_level=args.plots, outdir=args.outdir,
                        inlib=args.in_library, suffix=args.suffix)


def lincomb_spectrum(args):
    core.lincomb_spectrum(args.match_results, plot_level=args.plots,
                          inlib=args.in_library, outdir=args.outdir,
                          num_best=args.num_best, suffix=args.suffix)


def shift_spectrum(args):
    core.shift_spectrum(args.spectrum, indir=args.directory,
                        plot_level=args.plots, outdir=args.outdir,
                        suffix=args.suffix)


def main():
    args = sys.argv[1:]

    psr = ArgumentParser(description="SpecMatch-Emp: A tool for extracting " +
                         "fundamental stellar parameters from spectra",
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
    psr_sm.add_argument("-p", "--plots", action='count',
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
    psr_sm.add_argument("--n-lib-subset", type=int, default="",
                        help="Number of random stars to select from library. Useful for the purposes of debugging")

    psr_sm.set_defaults(func=specmatch_spectrum)

    psr_shift = subpsr.add_parser("shift")
    psr_shift.add_argument("spectrum", type=str, help="Path to spectrum file" +
                           " or cps observation id")
    psr_shift.add_argument("-d", "--directory", type=str,
                           default=os.path.join(SPECMATCHDIR, 'spectra'),
                           help="Directory to look in for spectra if an obs " +
                           "id was provided")
    psr_shift.add_argument("-o", "--outdir", type=str,
                           default=os.path.join(SPECMATCHDIR, 'shifted_spectra'),
                           help="Directory to store output files.")
    psr_shift.add_argument("-p", "--plots", action='count',
                           help="Generate diagnostic plots")
    psr_shift.add_argument("-s", "--suffix", type=str, default="_adj",
                        help="Suffix to append to results files")
    psr_shift.set_defaults(func=shift_spectrum)
    # psr_shift.add_argument("obs", type=str, help="cps id of target spectrum")
    # psr_shift.add_argument("-d", "--directory", type=str,
    #                        default=os.path.join(SPECMATCHDIR, 'spectra'),
    #                        help="Directory to look in for spectra")
    # psr_shift.add_argument("-n", "--name", type=str, default="",
    #                        help="Name to use as target ID")

    psr_match = subpsr.add_parser("match")
    psr_match.add_argument("spectrum", type=str, help="Path to spectrum file" +
                           " or cps observation id")
    psr_match.add_argument("-d", "--directory", type=str,
                           default=os.path.join(SPECMATCHDIR, 'spectra'),
                           help="Directory to look in for spectra if an obs " +
                           "id was provided")
    psr_match.add_argument("-p", "--plots", action='count',
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
    psr_lincomb.add_argument("-p", "--plots", action='count',
                             help="Generate diagnostic plots")
    psr_lincomb.add_argument("-i", "--in_library", type=str, default=False,
                             help="Check the library for a star with the " +
                             "same name and exclude it from the match.")
    psr_lincomb.add_argument("-o", "--outdir", type=str, default="./",
                             help="Directory to store output files.")
    psr_lincomb.add_argument("-n", "--num_best", type=int, default=5,
                             help="Number of best matches to use in linear " +
                             "combination step.")
    psr_lincomb.add_argument("-s", "--suffix", type=str, default="",
                             help="Suffix to append to results files")
    psr_lincomb.set_defaults(func=lincomb_spectrum)

    args = psr.parse_args(args=args)
    args.func(args)


if __name__ == '__main__':
    main()
