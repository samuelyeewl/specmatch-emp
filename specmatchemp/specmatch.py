"""
@filename specmatch.py

Class to carry out the specmatch
"""
import os
import lmfit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import astropy.io.fits as fits

from time import strftime

from specmatchemp import SPECMATCH_VERSION
from specmatchemp import spectrum
from specmatchemp.library import Library
from specmatchemp import shift
from specmatchemp import match
from specmatchemp import analysis
from specmatchemp import plots


class SpecMatch(object):
    """SpecMatch class to perform the SpecMatch routine.

    Begin by using the shift() method to shift the target spectrum onto the
    library wavelength scale. This uses a bootstrapping approach to find the
    best reference spectrum to use as a template for shifting.

    Next, the match() method runs a match against each library spectrum,
    computing a chi-squared value for each target-library pair.

    Finally, use the lincomb() method to interpolate between the best library
    matches to obtain the final derived parameters, stored in the results dict.

    Attributes:
        target (spectrum.Spectrum): Target spectrum object
        target_unshifted (spectrum.Spectrum): An unshifted spectrum
        shift_ref (spectrum.Spectrum): Reference used to shift the spectrum
        shift_results (dict): Dictionary containing results from shifting.
        regions (list of tuples): List of wavelength regions used for matching.
        match_results (pd.DataFrame): Parameter table including chi_squared
            and fit_params from single match routine.
        lincomb_matches (list of match.MatchLincomb): MatchLincomb objects for
            final matching to a linear combination of spectra
        lincomb_results (list of dicts):
            Results from each MatchLincomb fit.
        results (dict):
            Dictionary of final stellar parameters produced by SpecMatch.
            These are derived by taking the average of the detrended
            parameters from each wavelength region.
            Keys are elements of library.STAR_PROPS.

    Args:
        target (spectrum.Spectrum): Target spectrum and uncertainty
        lib (library.Library): Library object to match against
        wavlim (tuple, optional):
            Wavelength limits to perform matching on.
    """

    def __init__(self, target, lib, wavlim=None):
        if wavlim is None:
            self.wavlim = lib.wavlim
        else:
            self.wavlim = wavlim
        self.regions = self.wavlim

        # truncate target spectrum and library
        if isinstance(target, spectrum.HiresSpectrum):
            self.target = None
            self.target_unshifted = target
        else:
            self.target = target.cut(*self.wavlim)
            self.target_unshifted = target

        self.lib = lib.wav_cut(*self.wavlim)

        self._shifted = False
        self.shift_ref = None
        self.shift_results = {}
        self.regions = None
        self.match_results = pd.DataFrame()
        self.mt_lincomb = None
        self.results = {}

        return

    def shift(self):
        """Shift the target spectrum onto the library wavelength scale.

        Uses the Mgb triplet region (5120 - 5200 A), a section of spectrum
        containing much information, to determine which spectrum to use
        as the reference for shifting. It does this by initially shifting
        the target spectrum in that region to the allowed references and
        comparing the heights of the cross-correlation peaks. The reference
        with the highest median peak is used as the reference to shift the
        entire spectrum.

        Returns:
            spectrum.Spectrum: The shifted spectrum, which is also stored in
                self.target
        """
        shift_refs = self.lib.header['shift_refs']
        if 'nso' in shift_refs and self.lib.nso is None:
            raise ValueError("Error: Library did not contain NSO spectrum.")

        shift_specs = []
        for obs in shift_refs:
            if obs == 'nso':
                shift_specs.append(self.lib.nso)
            else:
                idx = self.lib.get_index(obs)
                shift_specs.append(self.lib.get_spectrum(idx))

        # Try to use the Mg triplet to determine which reference spectrum is
        # best.
        if isinstance(self.target_unshifted, spectrum.HiresSpectrum):
            # In the HIRES Echelle object, the Mgb triplet is in order 2
            ref_order = 2
            targ = self.target_unshifted
            targ_cut = spectrum.HiresSpectrum(targ.w[ref_order],
                                              targ.s[ref_order],
                                              targ.serr[ref_order])
        else:
            # If we already have a flattened spectrum, use the 5120-5200 A
            # region.
            targ_cut = self.target_unshifted.cut(5120, 5200)
            if len(targ_cut) == 0:
                # if 5120-5200 A region is missing, use the entire spectrum
                targ_cut = self.target_unshifted

        # Find the height of the correlation peak with each reference.
        median_peak = []
        for ref_spec in shift_specs:
            shift_data = {}
            shift.shift(targ_cut, ref_spec, store=shift_data)

            # get correlation peaks
            num_sects = shift_data['order_0/num_sections']
            peaks = []
            for sect in range(num_sects):
                xcorr = shift_data['order_0/sect_{0:d}/xcorr'.format(sect)]
                peaks.append(max(xcorr))

            median_peak.append(np.median(peaks))

        best_ref = shift_specs[np.argmax(median_peak)]

        # Now shift to the best reference
        self.target = shift.shift(self.target_unshifted, best_ref,
                                  store=self.shift_results)
        self.shift_ref = best_ref

        self._shifted = True
        return self.target

    def match(self, wavlim=None, wavstep=100, ignore=None):
        """Match the target against the library spectra.

        Performs a pairwise match between the target spectrum and every
        library spectra, fitting for vsini broadening and a spline to the
        continuum level. For each target-reference pair, a best chi-squared
        value is calculated to assess the match between the two spectra.

        The results are stored in
        :py:attr:`specmatchemp.SpecMatch.match_results`

        Args:
            wavlim (tuple or list of tuples, optional):
                Wavelength region(s) to use in matching process.
                Defaults to None - use the entire region overlapped by the
                target spectrum and library, split into sections specified
                by wavstep.
                If a tuple is given, it will be split into sections specified
                by wavstep.
                If a list of tuples is given, each tuple is a section
                to be used for matching.
            wavstep (float, optional):
                Length of wavelength regions to be used.
                Defaults to 100 Angstrom regions.
                If None, uses the entire region specified in wavlims
            ignore (int, optional): A library index to ignore. Useful for
                n-1 library test.
        """
        # Ensure spectrum has been shifted
        if self.target is None:
            print("Error: Target spectrum has not been shifted onto the " +
                  "library wavelength scale. Run SpecMatch.shift() first.")
            return

        # Create match results table
        self.match_results = self.lib.library_params.copy()

        # Get list of wavelength regions
        if isinstance(wavlim, list) or isinstance(wavlim, np.ndarray):
            regions = wavlim
        elif isinstance(wavlim, tuple) or wavlim is None:
            if wavlim is None:
                # If no wavlim is provided, use the library wavlim
                wavlim = self.wavlim

            if wavstep is None:
                # If no wavstep is provided, use the entire region given
                regions = [wavlim]
            else:
                # Split the region into sections
                startwl = np.arange(wavlim[0], wavlim[1], wavstep)
                regions = [(wl, wl + wavstep) for wl in startwl]
                # ensure final region doesn't exceed given bound
                regions[-1] = (regions[-1][0], wavlim[1])

        regions.sort()
        # ensure regions don't exceed either the spectrum or library bounds
        regions[0] = (max(regions[0][0], self.wavlim[0], self.target.w[0]),
                      regions[0][1])
        regions[-1] = (regions[-1][0], min(regions[-1][1], self.wavlim[1],
                          self.target.w[-1]))

        # save regions
        self.regions = regions

        for reg in regions:
            if len(regions) == 1:
                cs_col = 'chi_squared'
                fit_col = 'fit_params'
            else:
                cs_col = 'chi_squared_{0:.0f}'.format(reg[0])
                fit_col = 'fit_params_{0:.0f}'.format(reg[0])

            self.match_results.loc[:, cs_col] = np.nan
            self.match_results.loc[:, fit_col] = np.nan

            spec_targ = self.target.cut(*reg)

            for param_ref, spec_ref in self.lib:
                # ignore specified index
                if ignore is not None and param_ref.name == ignore:
                    continue

                # match
                mt = match.Match(spec_targ, spec_ref.cut(*reg))
                mt.best_fit()

                # store results
                ref_idx = param_ref.lib_index
                self.match_results.loc[ref_idx, cs_col] = mt.best_chisq
                self.match_results.loc[ref_idx, fit_col] = \
                    mt.best_params.dumps()

    def lincomb(self, num_best=5, regions='all'):
        """Interpolate between library spectra to get more accurate parameters.

        Takes the best `num_best` matches obtained in the match() step and
        creates linear combinations of their spectra. The respective
        coefficients for each spectrum will be used to take a weighted average
        of their parameters for use as the final derived parameters.

        Args:
            num_best (int, optional): Specify the number of best matches to be
                used to synthesize the linear combinations.
            regions (list of tuples, optional): Select wavelength regions to
                perform the matching procedure on.
                Defaults to 'all', which uses all the regions used in the
                previous match() step.
        """
        # Ensure matching process has already been run
        if self.match_results.empty:
            print("Error: Matching procedure has not yet been performed.\n" +
                  "Run SpecMatch.match() first.")
            return

        if regions == 'all':
            regions = self.regions
        elif not isinstance(regions, list):
            raise TypeError("regions should be a list of tuples")

        # Now perform lincomb match
        self.num_best = num_best
        self.ref_idxs = []
        self.vsini = []
        self.spec_refs = []
        for reg in regions:
            if len(regions) == 1:
                cs_col = 'chi_squared'
                fit_col = 'fit_params'
            else:
                cs_col = 'chi_squared_{0:.0f}'.format(reg[0])
                fit_col = 'fit_params_{0:.0f}'.format(reg[0])

            # get best matches
            self.match_results.sort_values(by=cs_col, inplace=True)
            ref_idxs = np.array(self.match_results.head(num_best).index)
            spec_refs = self.lib.get_spectrum(ref_idxs)
            spec_refs = [s.cut(*reg) for s in spec_refs]

            vsini = []
            for i in range(num_best):
                params = lmfit.Parameters()
                params.loads(self.match_results.iloc[i][fit_col])
                vsini.append(params['vsini'].value)
            vsini = np.array(vsini)

            mt_lincomb = match.MatchLincomb(self.target.cut(*reg),
                                            spec_refs, vsini)
            mt_lincomb.best_fit()
            coeffs = np.array(mt_lincomb.get_lincomb_coeffs())

            # obtain parameters
            lincomb_results = []
            lib_params = self.lib.library_params
            for p in Library.STAR_PROPS:
                lincomb_results[p] = analysis.lincomb_props(lib_params, p,
                                                            ref_idxs, coeffs)

            # save lincomb input and results
            self.ref_idxs.append(ref_idxs)
            self.spec_refs.append(spec_refs)
            self.lincomb_matches.append(mt_lincomb)
            self.lincomb_results.append(lincomb_results)

        # Average over all wavelength regions
        for p in Library.STAR_PROPS:
            self.results[p] = 0
            for i in range(len(regions)):
                self.results[p] += self.lincomb_results[i][p] / len(regions)

    def to_fits(self, outpath):
        """Saves the current state of the SpecMatch object to a fits file.

        Args:
            outpath (str): Path to store output file.
        """
        if not os.path.exists(os.path.dirname(outpath)):
            os.mkdir(os.path.dirname(outpath))

        hdulist = fits.HDUList()

        # Create Primary HDU.
        primhdu = fits.PrimaryHDU()
        primhdu.header['DATE'] = (strftime('%Y-%m-%d'),
                                  'File creation date (YYYY-MM-DD)')
        primhdu.header['VERSION'] = (SPECMATCH_VERSION,
                                     'SpecMatch-Emp version')
        primhdu.header['LIBDATE'] = (self.lib.header['date_created'],
                                     'Library creation date (YYYY-MM-DD)')
        primhdu.header.add_comment("SpecMatch-Emp Object", before="SIMPLE")

        hdulist.append(primhdu)

        # Save target spectrum
        if self.target is not None:
            targ_hdu = self.target.to_hdu()
            targ_hdu.name = 'TARGET'
            hdulist.append(targ_hdu)

        ################# Shift results #################
        # Save unshifted spectrum
        if self._shifted is True:
            # (Check if unshifted and shifted targets are different)
            unshifted_hdu = self.target_unshifted.to_hdu()
            unshifted_hdu.name = 'UNSHIFTED'
            hdulist.append(unshifted_hdu)

        # Save shift references
        if self.shift_ref is not None:
            shiftref_hdu = self.shift_ref.to_hdu()
            shiftref_hdu.name = 'SHIFTREF'
            hdulist.append(shiftref_hdu)

        # Save shift data
        if len(self.shift_results) > 0:
            shift_results = self.shift_results.copy()
            col_list = []
            num_orders = shift_results.pop('num_orders')
            col_list.append(fits.Column(name='num_orders', format='J',
                                        array=[num_orders]))

            num_sects = []
            for i in range(num_orders):
                num_sects.append(shift_results.pop('order_{0:d}/num_sections'
                                               .format(i)))
            col_list.append(fits.Column(name='num_sects', format='J',
                                        array=num_sects))
            n_sec = max(num_sects)

            for k in ['center_pix', 'lag', 'fit']:
                col_list.append(fits.Column(name=k,
                                format='{0:d}E'.format(n_sec),
                                array=shift_results.pop(k)))

            # Save individual fit data
            for k in shift_results.keys():
                col_list.append(fits.Column(name=k, format='D',
                                            array=shift_results[k]))

            shift_hdu = fits.BinTableHDU.from_columns(col_list)
            shift_hdu.name = 'SHIFTDATA'
            hdulist.append(shift_hdu)

        ################# Match Results #################
        if self.regions is not None:
            reg_hdu = fits.ImageHDU(name='REGIONS',
                                    data=np.array(self.regions))
            hdulist.append(reg_hdu)

        # Convert match results to record array
        if not self.match_results.empty:
            match_rec = self.match_results.to_records()
            dt = match_rec.dtype.descr
            for i in range(len(dt)):
                if dt[i][1] == "|O":
                    # max string length = 1000
                    dt[i] = (dt[i][0], 'S1000')
            match_rec = np.array(match_rec, dtype=dt)

            match_hdu = fits.BinTableHDU(name='MATCHRES',
                                         data=match_rec)
            hdulist.append(match_hdu)

        ################# Lincomb Results #################
        











    def plot_chi_squared_surface(self, num_best=None):
        """Plot the chi-squared surface from the pairwise matching procedure.

        Creates a three-column plot of the best chi-squared obtained with
        each library spectrum as a function of the library parameters.

        Args:
            num_best (optional [int]): Number of best spectra to highlight.
                (default: `self.num_best`)
        """
        if self.match_results.empty:
            return

        cs_col = 'chi_squared'
        if num_best is None:
            num_best = self.num_best

        ax1 = plt.subplot(131)
        plt.semilogy()
        plt.plot(self.match_results['Teff'], self.match_results[cs_col], '.')
        plt.plot(self.match_results['Teff'].head(num_best)\
            , self.match_results[cs_col].head(num_best), 'r.')
        plt.ylabel(r'$\chi^2$')
        plt.xlabel(r'$T_{eff}$ (K)')
        plt.xticks([3000,4000,5000,6000,7000])
        plots.reverse_x()

        ax2 = plt.subplot(132, sharey=ax1)
        plt.semilogy()
        plt.plot(self.match_results['radius'], self.match_results[cs_col], '.')
        plt.plot(self.match_results['radius'].head(num_best)\
            , self.match_results[cs_col].head(num_best), 'r.')
        ax2.set_xscale('log')
        plt.xlabel(r'$R\ (R_\odot)$')

        ax3 = plt.subplot(133, sharey=ax1)
        plt.semilogy()
        plt.plot(self.match_results['feh'], self.match_results[cs_col], '.')
        plt.plot(self.match_results['feh'].head(num_best)\
            , self.match_results[cs_col].head(num_best), 'r.')
        plt.xlabel(r'$[Fe/H]$ (dex)')

    def plot_best_matches(self, num_best=None):
        """Plots the reference, modified reference and residuals for each of the
        best matches.

        Args:
            num_best (optional [int]): Number of best spectra to highlight.
                (default: `self.num_best`)
        """
        if self.match_results is None:
            return

        if num_best is None:
            num_best = self.num_best

        fit_col = 'fit_params'

        shared_ax = None
        for i in range(num_best):
            if shared_ax is None:
                ax = plt.subplot(num_best, 1, i+1)
                shared_ax = ax
            else:
                ax = plt.subplot(num_best, 1, i+1, sharex=shared_ax)
            ref_idx = self.match_results.iloc[i]['lib_index.1']
            ref = self.lib.get_spectrum(ref_idx)
            params = lmfit.Parameters()
            params.loads(self.match_results.iloc[i][fit_col])
            mt = match.Match(self.target, ref)
            mt.load_params(params)
            mt.plot()
            if i != num_best-1:
                plt.setp(ax.get_xticklabels(), visible=False)
            ax.set_xlabel('')
            ax.set_ylabel('')

        # Label
        fig = plt.gcf()
        axes = fig.axes
        axes[-1].set_xlabel('Wavelength (Angstroms)', fontsize='large')
        fig.text(0.05, 0.5, 'Normalized Flux (Arbitrary Offset)', \
            ha='center', va='center', rotation='vertical', fontsize='large')
        plt.tight_layout(rect=[0.04,0,1,0.95])

    def plot_references(self, verbose=True, num_best=None):
        """Plots the location of the best references used in the linear
        combination step.

        Args:
            verbose (optional [bool]): Whether to annotate the points with
                the lincomb coefficients
            num_best (optional [int]): Number of best spectra to highlight.
                (default: `self.num_best`)
        """
        if self.match_results is None:
            return

        if num_best is None:
            num_best = self.num_best

        def _plot_ref_params(paramx, paramy, zoom=False):
            plt.plot(self.match_results[paramx], self.match_results[paramy], \
                '.', alpha=0.6, label='_nolegend_')
            plt.plot(self.match_results.head(num_best)[paramx], \
                self.match_results.head(num_best)[paramy], '^', color='forestgreen', label='Best Matches')
            if zoom:
                plots.set_tight_lims(self.match_results.head(num_best)[paramx], \
                    self.match_results.head(num_best)[paramy])
                if verbose: 
                    for i in range(self.num_best):
                        p = self.match_results.iloc[i]
                        plots.annotate_point(p[paramx], p[paramy], '{0:.3f}'.format(self.coeffs[i]))


        gs = gridspec.GridSpec(2,2)
        plt.subplot(gs[0])
        _plot_ref_params('Teff','radius')
        plt.plot(self.results['Teff'], self.results['radius'], 's', color='purple', label='Derived Parameters')
        ax = plt.gca()
        ax.set_yscale('log')
        plt.ylim((0.1, 16))
        plots.reverse_x()
        plt.xlabel(r'$T_{eff}$ (K)')
        plt.ylabel(r'$R\ (R_\odot)$')

        plt.subplot(gs[1])
        _plot_ref_params('Teff', 'radius', zoom=True)
        plt.plot(self.results['Teff'], self.results['radius'], 's', color='purple', label='Derived Parameters')
        ax = plt.gca()
        plots.reverse_x()
        plt.xlabel(r'$T_{eff}$ (K)')
        plt.ylabel(r'$R\ (R_\odot)$')

        plt.subplot(gs[2])
        _plot_ref_params('feh','radius')
        plt.plot(self.results['feh'], self.results['radius'], 's', color='purple', label='Derived Parameters')
        ax = plt.gca()
        ax.set_yscale('log')
        plt.ylim((0.1, 16))
        plt.xlabel(r'$[Fe/H]$ (dex)')
        plt.ylabel(r'$R\ (R_\odot)$')

        plt.subplot(gs[3])
        _plot_ref_params('feh', 'radius', zoom=True)
        plt.plot(self.results['feh'], self.results['radius'], 's', color='purple', label='Derived Parameters')
        ax = plt.gca()
        plt.xlabel(r'$[Fe/H]$ (dex)')
        plt.ylabel(r'$R\ (R_\odot)$')


    def plot_lincomb(self):
        """Plots the spectra used to synthesize the linear combination,
        the synthesized spectrum, and residuals.
        """
        if self.mt_lincomb is None:
            return

        self.mt_lincomb.plot()
        

