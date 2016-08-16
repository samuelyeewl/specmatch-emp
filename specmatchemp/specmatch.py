"""
@filename specmatch.py

Class to carry out the specmatch
"""
import lmfit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

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
        match_results (pd.DataFrame): Parameter table including chi_squared
            and fit_params from single match routine.
        mt_lincomb (match.MatchLincomb): MatchLincomb object for final
            matching to a linear combination of spectra
        results (dict): Dictionary of final stellar parameters produced
            by SpecMatch. Keys are elements of library.STAR_PROPS.

    Args:
        target (spectrum.Spectrum): Target spectrum and uncertainty
        lib (library.Library): Library object to match against
        wavlim (tuple, optional):
            Wavelength limits to perform matching on.
        num_best (optional [int]): Number of target spectra to use
            during linear combination matching stage.
    """

    def __init__(self, target, lib, wavlim=None, num_best=5):
        if wavlim is None:
            self.wavlim = lib.wavlim
        else:
            self.wavlim = wavlim

        # truncate target spectrum and library
        if isinstance(target, spectrum.HiresSpectrum):
            self.target = None
            self.target_unshifted = target
        else:
            self.target = target.cut(*wavlim)
            self.target_unshifted = target

        self.lib = lib.wav_cut(*wavlim)

        self.num_best = num_best
        self.shift_results = {}
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

        return self.target

    def match(self, ignore=None):
        """Match the target against the library spectra.

        Performs a pairwise match between the target spectrum and every
        library spectra, fitting for vsini broadening and a spline to the
        continuum level. For each target-reference pair, a best chi-squared
        value is calculated to assess the match between the two spectra.

        The results are stored in
        :py:attr:`specmatchemp.SpecMatch.match_results`

        Args:
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
        cs_col = 'chi_squared'
        fit_col = 'fit_params'
        self.match_results.loc[:, cs_col] = np.nan
        self.match_results.loc[:, fit_col] = np.nan

        for param_ref, spec_ref in self.lib:
            # ignore a particular index
            if ignore is not None and param_ref.name == ignore:
                continue

            # match
            mt = match.Match(self.target, spec_ref, opt='nelder')
            mt.best_fit()

            # store results
            ref_idx = param_ref.lib_index
            self.match_results.loc[ref_idx, cs_col] = mt.best_chisq
            self.match_results.loc[ref_idx, fit_col] = \
                mt.best_params.dumps()

        self.match_results.sort_values(by=cs_col, inplace=True)

        # if not lincomb:
        #     best_idx = self.match_results.iloc[0].name
        #     for p in Library.STAR_PROPS:
        #         self.results[p] = self.lib.library_params.loc[best_idx, p]
        #     return

        # # Now perform lincomb match
        # self.ref_idxs = np.array(self.match_results.head(self.num_best).index)
        # self.spec_refs = self.lib.get_spectrum(self.ref_idxs)
        # # get vsini
        # self.vsini = []
        # for i in range(self.num_best):
        #     params = lmfit.Parameters()
        #     params.loads(self.match_results.iloc[i][fit_col])
        #     self.vsini.append(params['vsini'].value)
        # self.vsini = np.array(self.vsini)

        # self.mt_lincomb = match.MatchLincomb(self.target, self.spec_refs,
        #                                      self.vsini)
        # self.mt_lincomb.best_fit()

        # # get derived values
        # self.coeffs = np.array(self.mt_lincomb.get_lincomb_coeffs())
        # self.get_derived_values()

        # return

    # get derived values
    def get_derived_values(self):
        for p in Library.STAR_PROPS:
            self.results[p] = analysis.lincomb_props(self.lib.library_params,
                                                p, self.ref_idxs, self.coeffs)

        return self.results

    def plot_chi_squared_surface(self, num_best=None):
        """Plot the chi-squared surface from the pairwise matching procedure.

        Creates a three-column plot of the best chi-squared obtained with
        each library spectrum as a function of the library parameters.

        Args:
            num_best (optional [int]): Number of best spectra to highlight.
                (default: `self.num_best`)
        """
        if self.match_results is None:
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
        

