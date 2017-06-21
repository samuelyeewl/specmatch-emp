"""
@filename core.py

SpecMatch-Emp core functions
"""
import os
import sys
from shutil import copy
import logging

import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

from specmatchemp import SPECMATCHDIR
from specmatchemp import SHIFT_REFERENCES
from specmatchemp import spectrum
from specmatchemp import shift
from specmatchemp import specmatch
from specmatchemp import library


def specmatch_spectrum(specpath, plot_level=0, inlib=False, outdir="./",
                       num_best=5, suffix="", wavlim='all', lib_subset=None,
                       name=None, n_lib_subset=None):
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

    if wavlim == 'all':
        target = spectrum.read_hires_fits(specpath)
    else:
        target = spectrum.read_hires_fits(specpath).cut(*wavlim)

    # Determine the name of the target
    if inlib:
        name = inlib
    elif name is None:
        name = os.path.basename(specpath)[:-5]

    if n_lib_subset is not None:
        lib = library.read_hdf(wavlim='none')
        lib_subset = lib.library_params.lib_index
        lib_subset = np.random.choice(
            lib_subset, size=n_lib_subset, replace=False
        )

    lib = library.read_hdf(wavlim=wavlim, lib_index_subset=lib_subset)
    sm = specmatch.SpecMatch(target, lib)
    sm.shift()

    if inlib:
        targ_idx = lib.get_index(inlib)
        targ_param, targ_spec = lib[targ_idx]
        sm.match(ignore=targ_idx, wavlim=wavlim)
    else:
        targ_param = None
        sm.match(wavlim=wavlim)

    sm.target.name = name  # attach target name

    sm.lincomb(num_best)

    # Print results
    print("SpecMatch Results for {0}".format(name))
    sm.results_to_txt(sys.stdout)

    outdir = os.path.join(outdir, name)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

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
        print("created {}".format(outpath))

    # Save full results
    outpath = os.path.join(outdir, name + suffix + '_sm.hdf')
    sm.to_hdf(outpath)
    print("created {}".format(outpath))

    # Create representative plots
    if plot_level is not None and plot_level > 0:
        plotspath = os.path.join(outdir, name + suffix + '_plots.pdf')
        with PdfPages(plotspath) as pdf:
            region = (5100, 5200)
            wavlim = (5160, 5190)
            order = 2

            plot_shifts(sm, pdf, order, wavlim)
            plot_match(sm, pdf, region, wavlim, targ_param)
            plot_lincomb(sm, pdf, region, wavlim, targ_param)
            print("created {}".format(plotspath))

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

        print("created {}".format(plotspath))

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
    fig, axL = plt.subplots(ncols=2, figsize=(10, 5))
    plt.sca(axL[0])
    sm.plot_xcorr(order, True)
    plt.title('{0} cross-correlation for order {1:d}'.format(name, order))
    plt.sca(axL[1])
    sm.plot_xcorr(order, True)
    meanshift = np.nanmean(sm.shift_data['fit'])
    dpix = 30
    plt.xlim(meanshift-dpix, meanshift+dpix)
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
    plt.tight_layout(rect=[0, 0, 1, 0.95])
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
    sm.plot_lincomb(region=region, wavlim=wavlim)
    plt.title('{0} Linear Combination results'.format(name))
    fig.set_tight_layout(True)
    pdf.savefig()
    plt.close()


def match_spectrum(specpath, indir="./", plot_level=0, inlib=False,
                   outdir="./", suffix=""):
    """Match a spectrum given its observation ID

    Args:
        specpath (str): Path to spectrum or its CPS observation id jXX.XXXX
        indir (str): Directory to look in for target spectrum
        plot_level (int, 0-2): Level of plotting to save
        inlib (str or False): String to search within library for to exclude
            from matching process
        outdir (str): Output file directory
        suffix (str): String to append to output file names

    Returns:
        specmatch.SpecMatch object
    """
    # Check if specpath is a path or an observation ID
    if os.path.exists(specpath):
        targ_path = specpath
        targid = os.path.splitext(os.path.basename(specpath))[0]
    else:
        targ_path = os.path.join(indir, 'r' + specpath + '_adj' + suffix +
                                 '.fits')
        if not os.path.exists(targ_path):
            raise ValueError(specpath + " does not exist!")
        targid = 'r' + specpath

    # Load shifted spectrum
    target = spectrum.read_fits(targ_path)

    lib = library.read_hdf()
    sm = specmatch.SpecMatch(target, lib)

    if inlib:
        name = inlib
        sm.target.name = inlib
        sm.target.attrs['obs'] = targid
        targ_idx = lib.get_index(inlib)
        targ_param, targ_spec = lib[targ_idx]
        sm.match(ignore=targ_idx)
    else:
        name = targid
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


def shift_spectrum(specpath, plot_level=0, indir=None, outdir="./",
                   suffix="_adj", mask=True, no_bootstrap=False,
                   flatten=False):
    """Shift a target spectrum given an observation code.

    Saves the shifted spectrum in a fits file.

    Args:
        specpath (str): Path to spectrum or its CPS observation id jXX.XXXX
        plot_level (int, 0-2): Level of plotting to save
        indir (str): Directory to look in for target spectrum
        outdir (str): Directory to store output files
        suffix (str): String to append to output file names
        mask (bool): Use a mask to remove telluric lines
        no_bootstrap (bool): Shift a spectrum without bootstrapping
        flatten (bool): If multiple chips are provided, flatten into a single
            spectrum file.
    Returns:
        shifted, unshifted, shift_data
    """
    # Check if specpath is a path or an observation ID
    if os.path.exists(specpath):
        targ_path = specpath
        targid = os.path.splitext(os.path.basename(specpath))[0]
    else:
        return _multishift_spectrum(specpath, plot_level, indir, outdir,
                                    suffix, mask, no_bootstrap, flatten)

    # if a different directory is provided, copy the file into specmatchemp
    # working directory
    specdir = os.path.join(SPECMATCHDIR, 'spectra')
    shiftedspecdir = os.path.join(SPECMATCHDIR, 'shifted_spectra')
    if indir != specdir and indir is not None:
        copy(targ_path, specdir)

    # load target and references
    if mask:
        maskfile = os.path.join(SPECMATCHDIR, 'hires_telluric_mask.csv')
    else:
        maskfile = None
    targ_spec = spectrum.read_hires_fits(targ_path, maskfile)

    if no_bootstrap:
        # Shift directly onto NSO spectrum
        ref_specs = [spectrum.read_fits(os.path.join(shiftedspecdir,
                     'nso_adj.fits'))]
        shift_data = {}
        print("Shifting directly against NSO spectrum.")
        shifted = shift.shift(targ_spec, ref_specs[0], store=shift_data)
        shift_data['shift_reference'] = 0
    else:
        # Shift spectrum onto boostrapped spectra
        ref_specs = [spectrum.read_fits(os.path.join(shiftedspecdir,
                     r[0] + '_adj.fits')) for r in SHIFT_REFERENCES]
        shift_data = {}
        shifted = shift.bootstrap_shift(targ_spec, ref_specs, store=shift_data)

    # Save shifted spectrum
    outpath = os.path.join(shiftedspecdir, targid + suffix + '.fits')
    shift.save_shift_to_fits(outpath, shifted, targ_spec, shift_data,
                             clobber=True)

    if outdir != shiftedspecdir:
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        copy(outpath, outdir)

    # Generate representative plots
    if plot_level == 1:
        plotfile = os.path.join(outdir, targid + "_shift_plots.pdf")
        print("Saving plots to " + plotfile)

        with PdfPages(plotfile) as pdf:
            # Get reference used
            shift_ref = ref_specs[shift_data['shift_reference']]
            # Plot single order
            plot_shift_data(targ_spec, shifted, shift_ref, shift_data, pdf, 2)
    # Generate individual plots for every order
    elif plot_level == 2:
        plotfile = os.path.join(outdir, targid + "_shift_plots.pdf")
        print("Saving plots to " + plotfile)
        with PdfPages(plotfile) as pdf:
            # Get reference used
            shift_ref = ref_specs[shift_data['shift_reference']]
            num_orders = shift_data['num_orders']
            for i in range(num_orders):
                plot_shift_data(targ_spec, shifted, shift_ref, shift_data,
                                pdf, i, singleorder=True)

    return shifted, targ_spec, shift_data


def _multishift_spectrum(cps_id, plot_level=0, indir=None, outdir="./",
                        suffix="_adj", mask=True, no_bootstrap=False,
                        flatten=False):
    """Helper function to shift multiple chips"""
    # If an observation id is given, search for all chips and shift each in
    # turn
    bj_path = os.path.join(indir, 'b' + cps_id + '.fits')
    if os.path.exists(bj_path):
        print("Shifting bj chip...")
        bj = shift_spectrum(bj_path, plot_level=plot_level, indir=indir,
                            outdir=outdir, suffix=suffix, mask=mask,
                            no_bootstrap=no_bootstrap, flatten=False)
    else:
        bj = None

    rj_path = os.path.join(indir, 'r' + cps_id + '.fits')
    if os.path.exists(rj_path):
        print("Shifting rj chip...")
        rj = shift_spectrum(rj_path, plot_level=plot_level, indir=indir,
                            outdir=outdir, suffix=suffix, mask=mask,
                            no_bootstrap=no_bootstrap, flatten=False)
    else:
        rj = None

    ij_path = os.path.join(indir, 'i' + cps_id + '.fits')
    if os.path.exists(ij_path):
        print("Shifting ij chip...")
        ij = shift_spectrum(ij_path, plot_level=plot_level, indir=indir,
                            outdir=outdir, suffix=suffix, mask=mask,
                            no_bootstrap=no_bootstrap, flatten=False)
    else:
        ij = None

    if bj is None and rj is None and ij is None:
        raise ValueError("No observations corresponding to " + cps_id +
                         "could be found in " + indir)

    specs_shifted = []
    specs_unshifted = []
    specs_sd = []
    chips = []
    if bj is not None:
        specs_shifted.append(bj[0])
        specs_unshifted.append(bj[1])
        specs_sd.append(bj[2])
        chips.append('bj')
    if rj is not None:
        specs_shifted.append(rj[0])
        specs_unshifted.append(rj[1])
        specs_sd.append(rj[2])
        chips.append('rj')
    if ij is not None:
        specs_shifted.append(ij[0])
        specs_unshifted.append(ij[1])
        specs_sd.append(ij[2])
        chips.append('ij')

    if flatten:
        # Combine all chips into a single file
        print("Flattening {0:d} spectra".format(len(specs_shifted)))
        shiftedspecdir = os.path.join(SPECMATCHDIR, 'shifted_spectra')
        nso = spectrum.read_fits(os.path.join(shiftedspecdir,
                                              'nso_adj.fits'))

        shifted = spectrum.Spectrum.combine_spectra(specs_shifted,
            nso.w, name=cps_id, prefixes=chips)

        unshifted = spectrum.HiresSpectrum.combine_spectra(specs_unshifted,
            name=cps_id, prefixes=chips)

        shift_data = {}
        for (i, sd) in enumerate(specs_sd):
            for k, v in sd.items():
                shift_data[chips[i] + '/' + k] = v

        # Save flattened spectrum
        outpath = os.path.join(shiftedspecdir, cps_id + suffix + '.fits')
        shifted.to_fits(outpath, clobber=True)

        if outdir != shiftedspecdir:
            if not os.path.exists(outdir):
                os.mkdir(outdir)
            copy(outpath, outdir)

        return shifted, unshifted, shift_data
    else:
        return specs_shifted, specs_unshifted, specs_sd
