"""
@filename specmatch.py

Class to carry out the specmatch
"""
import lmfit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from specmatchemp import library
from specmatchemp import match
from specmatchemp import analysis
from specmatchemp import plots


class SpecMatch(object):
    """SpecMatch class to perform the SpecMatch routine.

    Attributes:
        match_results (pd.DataFrame): Parameter table including chi_squared
            and fit_params from single match routine.
        mt_lincomb (match.MatchLincomb): MatchLincomb object for final
            matching to a linear combination of spectra
        results (dict): Dictionary of final stellar parameters produced
            by SpecMatch. Keys are elements of library.STAR_PROPS.

    Args:
        target (np.ndarray): Target spectrum and uncertainty
        lib (library.Library): Library object to match against
        wavlim (tuple): Wavelength limits to perform matching on.
        num_best (optional [int]): Number of target spectra to use
            during linear combination matching stage.
    """

    def __init__(self, target, lib, wavlim, num_best=5):
        self.wavlim = wavlim
        # truncate target spectrum and library
        wav = lib.wav
        wavidx, = np.where((wav >= wavlim[0]) & (wav <= wavlim[1]))
        idxmin = wavidx[0]
        idxmax = wavidx[-1]+1

        self.target = target.cut(*wavlim)
        self.lib = library.Library(lib.wav[idxmin:idxmax], lib.library_spectra[:,:,idxmin:idxmax],
                lib.library_params, lib.header, wavlim, lib.param_mask)

        self.num_best = num_best
        self.match_results = pd.DataFrame()
        self.mt_lincomb = None
        self.results = {}

    def match(self, match_results=None, lincomb=True):
        """Run the SpecMatch routine to obtain the parameters.

        First performs a pairwise match between the target spectrum and every
        library spectra, fitting for vsini broadening and a spline to the
        continuum level. For each target-reference pair, a best chi-squared
        value is calculated to assess the match between the two spectra.

        Using the best `num_best` spectra, synthesises linear combinations
        of these spectra with varying weights. The set of weights which
        minimize the chi-squared difference between this synthetic spectrum
        and the target is used as the set of weights in a weighted average
        of the library parameters to produce the derived parameters.

        Args:
            match_results (optional [pd.DataFrame]): Load in an existing 
                match results table.
            lincomb (optional [bool]): Whether to perform the linear 
                combination. If set to false, uses the best match as the
                derived results.
        """
        if match_results is not None:
            self.match_results = match_results
        else:
            # First, perform standard match
            self.match_results = self.lib.library_params.copy()
            cs_col = 'chi_squared'
            self.match_results.loc[:,cs_col] = np.nan
            fit_col = 'fit_params'
            self.match_results.loc[:,fit_col] = np.nan

            for param_ref, spec_ref in self.lib:
                # match
                mt = match.Match(self.target, spec_ref, opt='nelder')
                mt.best_fit()

                # store results
                ref_idx = param_ref.lib_index
                self.match_results.loc[ref_idx,cs_col] = mt.best_chisq
                self.match_results.loc[ref_idx,fit_col] = mt.best_params.dumps()

        self.match_results.sort_values(by=cs_col, inplace=True)

        if not lincomb:
            best_idx = self.match_results.iloc[0].name
            for p in library.STAR_PROPS:
                self.results[p] = self.lib.library_params.loc[best_idx, p]
            return

        # Now perform lincomb match
        self.ref_idxs = np.array(self.match_results.head(self.num_best).index)
        self.spec_refs = self.lib.get_spectrum(self.ref_idxs)
        # get vsini
        self.vsini = []
        for i in range(self.num_best):
            params = lmfit.Parameters()
            params.loads(self.match_results.iloc[i][fit_col])
            self.vsini.append(params['vsini'].value)
        self.vsini = np.array(self.vsini)

        self.mt_lincomb = match.MatchLincomb(self.target, self.spec_refs, self.vsini)
        self.mt_lincomb.best_fit()

        # get derived values
        self.coeffs = np.array(self.mt_lincomb.get_lincomb_coeffs())
        for p in library.STAR_PROPS:
            self.results[p] = analysis.lincomb_props(self.lib.library_params, p, self.ref_idxs, self.coeffs)


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
            ref_idx = self.match_results.iloc[i]['lib_index']
            ref = self.lib.get_spectrum(ref_idx)
            params = lmfit.Parameters()
            params.loads(self.match_results.iloc[i][fit_col])
            mt = match.Match(self.target, ref)
            mt.load_params(params)
            mt.plot()
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
                self.match_results.head(num_best)[paramy], '^', color='forestgreen')
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
        ax = plt.gca()
        ax.set_yscale('log')
        plt.ylim((0.1, 16))
        plots.reverse_x()
        plt.xlabel(r'$T_{eff}$ (K)')
        plt.ylabel(r'$R\ (R_\odot)$')

        plt.subplot(gs[1])
        _plot_ref_params('Teff', 'radius', zoom=True)
        ax = plt.gca()
        plots.reverse_x()
        plt.xlabel(r'$T_{eff}$ (K)')
        plt.ylabel(r'$R\ (R_\odot)$')

        plt.subplot(gs[2])
        _plot_ref_params('feh','radius')
        ax = plt.gca()
        ax.set_yscale('log')
        plt.ylim((0.1, 16))
        plt.xlabel(r'$[Fe/H]$ (dex)')
        plt.ylabel(r'$R\ (R_\odot)$')

        plt.subplot(gs[3])
        _plot_ref_params('feh', 'radius', zoom=True)
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
        

