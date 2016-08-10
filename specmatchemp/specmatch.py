"""
@filename specmatch.py

Class to carry out the specmatch
"""
import lmfit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from specmatchemp import library
from specmatchemp import match
from specmatchemp import analysis
from specmatchemp.plotting import plots

class SpecMatch(object):
    """SpecMatch class to perform the SpecMatch routine.

    Attributes:
        match_results (pd.DataFrame): Parameter table including chi_squared
            and fit_params from Match routine.

    Args:
        target (np.ndarray): Target spectrum and uncertainty
        lib (library.Library): Library object to match against
        wavlim (tuple): Wavelength limits to perform matching on.
    """

    def __init__(self, target, lib, wavlim, num_best=5):
        self.wavlim = wavlim
        # truncate target spectrum and library
        wav = lib.wav
        wavidx, = np.where((wav >= wavlim[0]) & (wav <= wavlim[1]))
        idxmin = wavidx[0]
        idxmax = wavidx[-1]+1

        self.target = target[:,idxmin:idxmax]
        self.lib = library.Library(lib.wav[idxmin:idxmax], lib.library_spectra[:,:,idxmin:idxmax]\
            , lib.library_params, lib.header, wavlim, lib.param_mask)

        self.num_best = num_best
        self.match_results = pd.DataFrame()
        self.mt_lincomb = None
        self.results = {}

    def match(self):
        # First, perform standard match
        self.match_results = self.lib.library_params.copy()
        cs_col = 'chi_squared'
        self.match_results.loc[:,cs_col] = np.nan
        fit_col = 'fit_params'
        self.match_results.loc[:,fit_col] = np.nan

        for param_ref, spec_ref in self.lib:
            # match
            mt = match.Match(self.lib.wav, self.target, spec_ref, opt='nelder')
            mt.best_fit()

            # store results
            ref_idx = param_ref.lib_index
            self.match_results.loc[ref_idx,cs_col] = mt.best_chisq
            self.match_results.loc[ref_idx,fit_col] = mt.best_params.dumps()

        # Now perform lincomb match
        self.match_results.sort_values(by=cs_col, inplace=True)
        ref_idxs = np.array(self.match_results.head(self.num_best).index)
        spec_refs = self.lib.library_spectra[ref_idxs]
        # get vsini
        vsini = []
        for i in range(self.num_best):
            params = lmfit.Parameters()
            params.loads(self.match_results.iloc[i][fit_col])
            vsini.append(params['vsini'].value)
        vsini = np.array(vsini)

        self.mt_lincomb = match.MatchLincomb(self.lib.wav, self.target, spec_refs, vsini)
        self.mt_lincomb.best_fit()

        # get derived values
        coeffs = np.array(match.get_lincomb_coeffs(self.mt_lincomb.best_params))
        for p in library.STAR_PROPS:
            self.results[p] = analysis.lincomb_props(self.lib.library_params, p, ref_idxs, coeffs)

    def plot_chi_squared_surface(self):
        plt.subplot(131)
        plt.semilogy()
        plots.plot_param_chi_squared(self.match_results, 'Teff')
        plt.ylabel(r'$\chi^2$')
        plt.xlabel(r'$T_{eff}$ (K)')
        plt.xticks([3000,4000,5000,6000,7000])
        plots.reverse_x()
        plt.subplot(132)
        plt.semilogy()
        plots.plot_param_chi_squared(self.match_results, 'radius')
        ax = plt.gca()
        ax.set_xscale('log')
        plt.xlabel(r'$R\ (R_\odot)$')
        plt.subplot(133)
        plt.semilogy()
        plots.plot_param_chi_squared(self.match_results, 'feh')
        plt.xlabel(r'$[Fe/H]$ (dex)')

    # def plot_references(self):
    #     plt.plot(self.)

