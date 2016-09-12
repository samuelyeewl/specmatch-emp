"""
@filename core.py

SpecMatch-Emp core functions
"""

import os
from shutil import copy

import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

from specmatchemp import SPECMATCHDIR
from specmatchemp import SHIFT_REFERENCES
from specmatchemp import spectrum
from specmatchemp import shift
from specmatchemp import plots
from specmatchemp import specmatch
from specmatchemp import library


def specmatch_spectrum(specpath, plot_level=0, inlib=False, outdir="./",
                       num_best=5, suffix=""):
    """Perform the specmatch on a given spectrum

    Args:
        specpath (str): Path to target spectrum
        plot_level (int, 0-2): Level of plots
            0 - No plots saved, 1 - Representative plots, 2 - All plots
        inlib (str or False): String to search within library for to exclude
            from matching process
        outdir (str): Output file directory
        num_best (int): Number of best matches to use at lincomb stage
        suffix (str): String to append to output file names

    Returns:
        specmatch.SpecMatch object
    """
    if not os.path.exists(specpath):
        raise ValueError(specpath + " does not exist!")

    target = spectrum.read_hires_fits(specpath)
    lib = library.read_hdf()
    sm = specmatch.SpecMatch(target, lib)
    sm.shift()

    if inlib:
        sm.target.name = inlib
        name = inlib
        targ_idx = lib.get_index(inlib)
        targ_param, targ_spec = lib[targ_idx]
        sm.match(ignore=targ_idx)
    else:
        name = os.path.basename(specpath)[:-4]
        sm.target.name = name
        targ_param = None
        sm.match()

    sm.lincomb(num_best)

    # Print results
    print("SpecMatch Results for {0}".format(name))
    for p in library.Library.STAR_PROPS:
        print("{0}: {1:.2f}".format(p, sm.results[p]))

    outdir = os.path.join(outdir, name)
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # Save final results
    outpath = os.path.join(outdir, name + suffix + '_results.txt')
    with open(outpath, 'w') as f:
        if inlib:
            f.write('Library Parameters\n')
            f.write('------------------\n')
            f.write('Teff: {0:.0f} +/- {1:.0f} K\n'.format(
                targ_param['Teff'], targ_param['u_Teff']))
            f.write('Radius: {0:.3f} +/- {1:.3f} Rsun\n'.format(
                targ_param['radius'], targ_param['u_radius']))
            f.write('[Fe/H]: {0:.2f} +/- {1:.2f} dex\n'.format(
                targ_param['feh'], targ_param['u_feh']))
        f.write('\n')
        sm.results_to_txt(f, verbose=True)

    # Save full results
    outpath = os.path.join(outdir, name + suffix + '_sm.hdf')
    sm.to_hdf(outpath)

    # Create representative plots
    if plot_level > 0:
        plotspath = os.path.join(outdir, name + suffix + '_plots.pdf')
        with PdfPages(plotspath) as pdf:
            region = (5100, 5200)
            wavlim = (5160, 5190)
            order = 2

            plot_shifts(sm, pdf, order, wavlim)
            plot_match(sm, pdf, region, wavlim, targ_param)
            plot_lincomb(sm, pdf, region, wavlim, targ_param)

    # Create full plots
    if plot_level == 2:
        shiftplotspath = os.path.join(outdir, name + suffix +
                                      '_shift_plots.pdf')
        with PdfPages(shiftplotspath) as pdf:
            for order in range(np.shape(sm.target_unshifted.w)[0]):
                plot_shifts(sm, pdf, order, wavlim='all')

        matchplotspath = os.path.join(outdir, name + suffix +
                                      '_match_plots.pdf')
        with PdfPages(matchplotspath) as pdf:
            for reg in sm.regions:
                plot_match(sm, pdf, reg, wavlim='all', targ_param=targ_param)

        lincombplotspath = os.path.join(outdir, name + suffix +
                                        '_lincomb_plots.pdf')
        with PdfPages(lincombplotspath) as pdf:
            for reg in sm.lincomb_regions:
                plot_lincomb(sm, pdf, reg, wavlim='all', targ_param=targ_param)

    return sm


def plot_shifts(sm, pdf, order, wavlim='all', singleorder=False):
    """Create shift plots

    Args:
        sm (specmatch.SpecMatch): SpecMatch object to plot
        pdf: Pdf file object
        order (int): HIRES order to plot
        wavlim: A specific wavelength range to plot spectrum
        singleorder (bool): Whether to plot lags as a single order
    """
    name = sm.target.name
    if wavlim == 'all':
        min_w = sm.target_unshifted.w[order][0]
        max_w = sm.target_unshifted.w[order][-1]
        wavlim = (min_w, max_w)

    # Shifted spectrum
    fig = plt.figure(figsize=(10, 5))
    sm.plot_shifted_spectrum(wavlim=wavlim)
    plt.title('Shift results for star {0}'.format(name))
    fig.set_tight_layout(True)
    pdf.savefig()
    plt.close()

    # Cross-correlation
    fig = plt.figure(figsize=(10, 5))
    sm.plot_xcorr(order, True)
    plt.title('{0} cross-correlation for order {1:d}'
              .format(name, order))
    fig.set_tight_layout(True)
    pdf.savefig()
    plt.close()

    # Lags
    fig = plt.figure(figsize=(8, 6))
    if singleorder:
        sm.plot_shift_lags(order)
        plt.title('{0} lags for order {1:d}'.format(name, order))
    else:
        sm.plot_shift_lags()
        plt.title('{0} lags'.format(name))
    fig.set_tight_layout(True)
    pdf.savefig()
    plt.close()


def plot_shift_data(target_unshifted, target_shifted, reference, shift_data,
                    pdf, order, wavlim='all', singleorder=False):
    """Create shift plots from shift data

    Args:
        target_unshifted (HiresSpectrum): Unshifted target spectrum
        target_shifted (Spectrum): Shifted target spectrum
        reference (Spectrum): Reference spectrum
        shift_data (dict): Shift data object
        pdf: Pdf file object
        order (int): HIRES order to plot
        wavlim: A specific wavelength range to plot spectrum
    """
    # Create temp specmatch object
    sm = specmatch.SpecMatch(target_unshifted)
    sm.target = target_shifted
    sm.shift_ref = reference
    sm.shift_data = shift_data

    plot_shifts(sm, pdf, order, wavlim, singleorder)


def plot_match(sm, pdf, region=0, wavlim='all', targ_param=None):
    """Create match plots

    Args:
        sm (specmatch.SpecMatch): SpecMatch object to plot
        pdf: Pdf file object
        region (int or tuple): Match region to plot
        wavlim: A specific wavelength range to plot spectrum
        targ_param: Target parameters
    """
    name = sm.target.name

    # Chi-squared surfaces
    fig = plt.figure(figsize=(12, 8))
    sm.plot_chi_squared_surface()
    if targ_param is not None:
        # Plot target parameters if available
        axes = fig.axes
        axes[0].axvline(targ_param['Teff'], color='k')
        axes[1].axvline(targ_param['radius'], color='k')
        axes[2].axvline(targ_param['feh'], color='k')
    plt.suptitle('{0} chi-squared surface'.format(name))
    fig.set_tight_layout(True)
    pdf.savefig()
    plt.close()

    # Best match spectrum
    fig = plt.figure(figsize=(10, 4))
    sm.plot_best_match_spectra(region=region, wavlim=wavlim, num_best=1)
    plt.suptitle('{0} best matching spectrum'.format(name))
    plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])
    pdf.savefig()
    plt.close()


def plot_lincomb(sm, pdf, region=0, wavlim='all', targ_param=None):
    """Create lincomb plots

    Args:
        sm (specmatch.SpecMatch): SpecMatch object to plot
        pdf: Pdf file object
        region (int or tuple): Match region to plot
        wavlim: A specific wavelength range to plot spectrum
        targ_param: Target parameters
    """
    name = sm.target.name

    # Reference locations
    fig = plt.figure(figsize=(10, 8))
    sm.plot_references(region=region, verbose=True)
    axes = fig.axes
    axes[0].legend(numpoints=1, fontsize='small', loc='best')
    if targ_param is not None:
        # Plot target parameters if available
        axes[0].plot(targ_param['Teff'], targ_param['radius'], '*',
                     ms=15, color='red', label='Target')
        axes[1].plot(targ_param['Teff'], targ_param['radius'], '*',
                     ms=15, color='red')
        axes[2].plot(targ_param['feh'], targ_param['radius'], '*',
                     ms=15, color='red')
        axes[3].plot(targ_param['feh'], targ_param['radius'], '*',
                     ms=15, color='red')
    plt.suptitle('{0} references used in linear combination'
                 .format(name))
    plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])
    pdf.savefig()
    plt.close()

    # Lincomb results
    fig = plt.figure(figsize=(12, 6))
    sm.plot_lincomb(region=(5100, 5200), wavlim=(5160, 5190))
    plt.title('{0} Linear Combination results'.format(name))
    fig.set_tight_layout(True)
    pdf.savefig()
    plt.close()


def match_spectrum(obs, indir="./", plot_level=0, inlib=False, outdir="./",
                   suffix=""):
    """Match a spectrum given its observation ID

    Args:
        obs (str): CPS id of target spectrum
        indir (str): Directory to look in for target spectrum
        plot_level (int, 0-2): Level of plotting to save
        inlib (str or False): String to search within library for to exclude
            from matching process
        outdir (str): Output file directory
        suffix (str): String to append to output file names

    Returns:
        specmatch.SpecMatch object
    """
    inpath = os.path.join(indir, 'r' + obs + '_adj' + suffix + '.fits')
    if not os.path.exists(inpath):
        raise ValueError(inpath + " does not exist!")

    # Load shifted spectrum
    target = spectrum.read_fits(inpath)

    lib = library.read_hdf()
    sm = specmatch.SpecMatch(target, lib)

    if inlib:
        name = inlib
        sm.target.name = inlib
        sm.target.attrs['obs'] = 'r' + obs
        targ_idx = lib.get_index(inlib)
        targ_param, targ_spec = lib[targ_idx]
        sm.match(ignore=targ_idx)
    else:
        name = 'r' + obs
        targ_param = None
        sm.match()

    # Save results
    outdir = os.path.join(outdir, name)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    outpath = os.path.join(outdir, name + suffix + '_match.csv')
    sm.match_results.to_csv(outpath)

    # Save SpecMatch object
    outpath = os.path.join(outdir, name + suffix + '_sm.hdf')
    sm.to_hdf(outpath)

    # Generate representative plots
    if plot_level == 1:
        plotspath = os.path.join(outdir, name + suffix + '_match_plots.pdf')
        with PdfPages(plotspath) as pdf:
            region = (5100, 5200)
            wavlim = (5160, 5190)
            plot_match(sm, pdf, region, wavlim, targ_param)

    # Generate full plots
    if plot_level == 2:
        plotspath = os.path.join(outdir, name + suffix + '_match_plots.pdf')
        with PdfPages(plotspath) as pdf:
            for reg in sm.regions:
                plot_match(sm, pdf, reg, wavlim='all', targ_param=targ_param)

    return sm


def lincomb_spectrum(respath, plot_level=0, inlib=False, outdir="./",
                     num_best=5, suffix=""):
    """Match a spectrum using the linear combination approach.
    Can only be used to resume an existing SpecMatch object.

    Args:
        respath (str): Path to existing SpecMatch.hdf file
        plot_level (int, 0-2): Level of plotting to save
        inlib (str or False): String to search within library for to exclude
            from matching process
        outdir (str): Output file directory
        num_best (int): Number of best matches to use at lincomb stage
        suffix (str): String to append to output file names

    Returns:
        specmatch.SpecMatch object
    """
    lib = library.read_hdf()
    sm = specmatch.SpecMatch.read_hdf(respath, lib)
    name = sm.target.name

    if inlib:
        targ_idx = lib.get_index(inlib)
        targ_param, targ_spec = lib[targ_idx]
    else:
        targ_param = None

    sm.lincomb(num_best=num_best)

    # Print results
    print("SpecMatch Results for {0}".format(name))
    for p in library.Library.STAR_PROPS:
        print("{0}: {1:.2f}".format(p, sm.results[p]))

    outdir = os.path.join(outdir, name)
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # Save final results
    outpath = os.path.join(outdir, name + suffix + '_results.txt')
    with open(outpath, 'w') as f:
        if inlib:
            f.write('Library Parameters\n')
            f.write('------------------\n')
            f.write('Teff: {0:.0f} +/- {1:.0f} K\n'.format(
                targ_param['Teff'], targ_param['u_Teff']))
            f.write('Radius: {0:.3f} +/- {1:.3f} Rsun\n'.format(
                targ_param['radius'], targ_param['u_radius']))
            f.write('[Fe/H]: {0:.2f} +/- {1:.2f} dex\n'.format(
                targ_param['feh'], targ_param['u_feh']))
        f.write('\n')
        sm.results_to_txt(f, verbose=True)

    # Save full results
    outpath = os.path.join(outdir, name + suffix + '_lincomb_sm.hdf')
    sm.to_hdf(outpath)

    # Plot results
    if plot_level == 1:
        plotspath = os.path.join(outdir, name + suffix + '_lincomb_plots.pdf')
        with PdfPages(plotspath) as pdf:
            region = (5100, 5200)
            wavlim = (5160, 5190)
            plot_lincomb(sm, pdf, region, wavlim, targ_param)

    if plot_level == 2:
        plotspath = os.path.join(outdir, name + suffix + '_lincomb_plots.pdf')
        with PdfPages(plotspath) as pdf:
            for reg in sm.lincomb_regions:
                plot_lincomb(sm, pdf, reg, wavlim='all', targ_param=targ_param)
    return sm


def shift_spectrum(obs, indir="./", plot_level=0, outdir="./",
                   name="", suffix="", mask=True):
    """Shift a target spectrum given an observation code.

    Saves the shifted spectrum in a fits file.

    Args:
        obs (str): CPS id of target spectrum
        indir (str): Directory to look in for target spectrum
        plot_level (int, 0-2): Level of plotting to save
        name (str): Name to use as target ID.
        suffix (str): String to append to output file names
    Returns:
        shifted, unshifted, shift_data
    """
    # if a different directory is provided, copy the file into specmatchemp
    # working directory
    specdir = os.path.join(SPECMATCHDIR, 'spectra')
    shiftedspecdir = os.path.join(SPECMATCHDIR, 'shifted_spectra')

    if indir is not specdir:
        copy(os.path.join(indir, 'r' + obs + '.fits'), specdir)

    # load target and references
    targ_path = os.path.join(specdir, 'r' + obs + '.fits')
    if mask:
        maskfile = os.path.join(SPECMATCHDIR, 'hires_telluric_mask.csv')
    else:
        maskfile = None
    targ_spec = spectrum.read_hires_fits(targ_path, maskfile)

    if len(name) > 0:
        targ_spec.name = name
        targid = name
    else:
        targid = obs

    ref_specs = [spectrum.read_fits(os.path.join(shiftedspecdir,
                 r[0] + '_adj.fits')) for r in SHIFT_REFERENCES]

    # Shift spectrum
    shift_data = {}
    shifted = shift.bootstrap_shift(targ_spec, ref_specs, store=shift_data)
    # Save shifted spectrum
    outpath = os.path.join(shiftedspecdir, 'r' + obs + '_adj' +
                           suffix + '.fits')
    shift.save_shift_to_fits(outpath, shifted, targ_spec, shift_data,
                             clobber=True)
    if outdir is not shiftedspecdir:
        outdir = os.path.join(outdir, name)
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        copy(outpath, outdir)

    # Generate representative plots
    if plot_level == 1:
        plotfile = os.path.join(outdir, targid + "_shift_plots.pdf")
        wavlim = (5160, 5190)

        with PdfPages(plotfile) as pdf:
            # Get reference used
            shift_ref = ref_specs[shift_data['shift_reference']]
            # Plot single order
            plot_shift_data(targ_spec, shifted, shift_ref, shift_data, pdf, 2)
    # Generate individual plots for every order
    elif plot_level == 2:
        plotfile = os.path.join(outdir, targid + "_shift_plots.pdf")
        with PdfPages(plotfile) as pdf:
            # Get reference used
            shift_ref = ref_specs[shift_data['shift_reference']]
            num_orders = shift_data['num_orders']
            for i in range(num_orders):
                plot_shift_data(targ_spec, shifted, shift_ref, shift_data,
                                pdf, i, singleorder=True)


    return shifted, targ_spec, shift_data
