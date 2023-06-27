"""
@filename specmatch.py

Class to carry out the specmatch
"""
from six import string_types

import os
import csv
import lmfit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import h5py
import astropy.io.fits as fits

from time import strftime

from specmatchemp import SPECMATCHDIR
from specmatchemp import SPECMATCH_VERSION
from specmatchemp import SHIFT_REFERENCES
from specmatchemp.io import h5plus
from specmatchemp import spectrum
from specmatchemp.spectrum import Spectrum
from specmatchemp.spectrum import HiresSpectrum
from specmatchemp.library import Library
from specmatchemp import shift
from specmatchemp import match
from specmatchemp import analysis
from specmatchemp import plots
from specmatchemp import detrend

WAVLIM_DEFAULT = (5000, 5900)
WAVSTEP_DEFAULT = 100


class SpecMatch(object):
    """SpecMatch class to perform the SpecMatch routine.

    Begin by using the shift() method to shift the target spectrum onto the
    library wavelength scale. This uses a bootstrapping approach to find the
    best reference spectrum to use as a template for shifting.

    Next, the match() method runs a match against each library spectrum,
    computing a chi-squared value for each target-library pair.

    Finally, use the lincomb() method to interpolate between the best library
    matches to obtain the final derived parameters, stored in the results dict.

    Args:
        target (spectrum.Spectrum): Target spectrum and uncertainty
        lib (library.Library): Library object to match against
        wavlim (tuple, optional):
            Wavelength limits to perform matching on.

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
    """

    def __init__(self, target, lib=None, wavlim=WAVLIM_DEFAULT):
        if wavlim is None:
            self.wavlim = lib.wavlim
        elif lib is None:
            self.wavlim = wavlim
        else:
            self.wavlim = (max(wavlim[0], lib.wavlim[0]),
                           min(wavlim[1], lib.wavlim[1]))
        self.regions = self.wavlim

        # truncate target spectrum and library
        if isinstance(target, HiresSpectrum):
            self.target = None
            self.target_unshifted = target
        else:
            self.target = target.cut(*self.wavlim)
            self.target_unshifted = target

        if lib is None:
            self.lib = Library()
        else:
            self.lib = lib

        self._shifted = False
        self.num_best = 5
        self.shift_ref = None
        self.shift_data = {}
        self.regions = None
        self.lincomb_regions = None
        self.match_results = pd.DataFrame()
        self.lincomb_matches = []
        self.coeffs = None
        self.lincomb_results = []
        self.results_nodetrend = {}
        self.results = {}
        self.u_table = None

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
        print("Shifting spectrum")
        shift_refs = SHIFT_REFERENCES

        # Obtain spectra
        shift_specs = []
        for ref in shift_refs:
            obs = ref[0]
            if obs == 'nso':
                if self.lib.nso is None:
                    raise ValueError("Error: Library did not contain " +
                                     "NSO spectrum.")
                shift_specs.append(self.lib.nso)
            else:
                idx = self.lib.get_index('r'+obs)
                if idx is not None:
                    shift_specs.append(self.lib.get_spectrum(idx))

        # Use the bootstrap shift approach
        self.target = shift.bootstrap_shift(self.target_unshifted, shift_specs,
                                            store=self.shift_data)

        self.shift_ref = shift_specs[self.shift_data['shift_reference']]

        self._shifted = True
        return self.target

    def match(self, wavlim=WAVLIM_DEFAULT, wavstep=WAVSTEP_DEFAULT,
              ignore=None):
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

        print("Matching spectrum")
        # Create match results table
        self.match_results = self.lib.library_params.copy()

        # Get list of wavelength regions
        if isinstance(wavlim, list) or isinstance(wavlim, np.ndarray):
            regions = wavlim
        elif isinstance(wavlim, tuple) or wavlim is None or wavlim == 'all':
            if wavlim is None or wavlim == 'all':
                # If no wavlim is provided, use the default wavlim
                wavlim = WAVLIM_DEFAULT

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
        regions[0] = (max(regions[0][0], self.wavlim[0]), regions[0][1])
        regions[-1] = (regions[-1][0],
                       min(regions[-1][1], self.wavlim[1]))

        # save regions
        self.regions = regions

        for reg in regions:
            print("Matching region {0}".format(reg))
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
        print("Creating linear combinations")
        # Ensure matching process has already been run
        if self.match_results.empty:
            print("Error: Matching procedure has not yet been performed.\n" +
                  "Run SpecMatch.match() first.")
            return

        if regions == 'all':
            lincomb_regions = self.regions
        elif isinstance(regions, tuple):
            lincomb_regions = [regions]
        elif isinstance(regions, list):
            lincomb_regions = regions
        else:
            raise TypeError("regions should be a list of tuples")

        # Now perform lincomb match
        self.num_best = num_best
        self.ref_idxs = []
        self.coeffs = []
        self.lincomb_matches = []
        self.lincomb_results = []
        self.lincomb_regions = lincomb_regions
        for reg in lincomb_regions:
            print("Linear combinations in region {0}".format(reg))
            if len(self.regions) == 1:
                cs_col = 'chi_squared'
                fit_col = 'fit_params'
            else:
                cs_col = 'chi_squared_{0:.0f}'.format(reg[0])
                fit_col = 'fit_params_{0:.0f}'.format(reg[0])

            # get best matches
            self.match_results.sort_values(by=cs_col, inplace=True)
            ref_idxs = np.array(self.match_results.head(num_best).index)
            ref_chisq = np.array(self.match_results.head(num_best)[cs_col])
            spec_refs = self.lib.get_spectrum(ref_idxs)
            spec_refs = [s.cut(*reg) for s in spec_refs]

            vsini = []
            for i in range(num_best):
                params = lmfit.Parameters()
                params.loads(self.match_results.iloc[i][fit_col])
                vsini.append(params['vsini'].value)
            vsini = np.array(vsini)

            mt_lincomb = match.MatchLincomb(self.target.cut(*reg),
                                            spec_refs, vsini,
                                            ref_chisq=ref_chisq)
            mt_lincomb.best_fit()
            coeffs = np.array(mt_lincomb.get_lincomb_coeffs())

            # obtain parameters
            lincomb_results = {}
            lib_params = self.lib.library_params
            for p in Library.STAR_PROPS:
                lincomb_results[p] = analysis.lincomb_props(lib_params, p,
                                                            ref_idxs, coeffs)

            # save lincomb input and results
            self.ref_idxs.append(ref_idxs)
            self.lincomb_matches.append(mt_lincomb)
            self.lincomb_results.append(lincomb_results)
            self.coeffs.append(coeffs)

        # Read in uncertainties
        self._read_uncertainties()        
            
        # Average over all wavelength regions
        for p in Library.STAR_PROPS:
            self.results_nodetrend[p] = 0
            for i in range(len(lincomb_regions)):
                self.results_nodetrend[p] += (self.lincomb_results[i][p] /
                                              len(lincomb_regions))
            # Add uncertainties
            self.results_nodetrend['u_'+p] = self._get_uncertainty(self.results_nodetrend[p], p)

        # Detrend parameters
        d = detrend.Detrend()
        for p in Library.STAR_PROPS:
            self.results[p] = d.detrend(self.results_nodetrend[p], p)
            # TODO: Add uncertainties
            self.results['u_'+p] = self._get_uncertainty(self.results[p], p)

    def _read_uncertainties(self, filename=""):
        # Read in the uncertainties
        if len(filename) == 0:
            filename = os.path.join(SPECMATCHDIR, 'uncertainties.csv')

        self.u_table = {}

        with open(filename, mode='r') as csvfile:
            dialect = csv.Sniffer().sniff(csvfile.read(1024))
            csvfile.seek(0)

            reader = csv.reader(csvfile, dialect)
            # Skip header row
            next(reader)

            for row in reader:
                param = row[0]
                if param in self.u_table:
                    self.u_table[param].append((float(row[1]),
                        float(row[2]), float(row[3])))
                else:
                    self.u_table[param] = [(float(row[1]),
                        float(row[2]), float(row[3]))]

        for p in self.u_table:
            self.u_table[param].sort()

    def _get_uncertainty(self, value, param):
        if self.u_table is None or param not in self.u_table:
            return 0.0
        # Find appropriate interval
        for row in self.u_table[param]:
            if value >= row[1] and value < row[2]:
                return row[0]
        return 0.0

    @classmethod
    def read_hdf(cls, infile, lib):
        """Reads a SpecMatch object from and hdf file.

        Args:
            infile (str or h5 file): Input path or file handle.
            lib (library.Library): Library used to create SpecMatch object.
        """
        is_path = False
        if isinstance(infile, string_types):
            infile = h5py.File(infile, 'r')
            is_path = True

        # Read target spectrum
        target = spectrum.read_hdf(infile['target'])

        sm = cls(target, lib)

        # ------------------------ Shift results ------------------------
        if 'unshifted' in infile:
            sm.target_unshifted = spectrum.read_hdf(infile['unshifted'])

        if 'shift_ref' in infile:
            sm.shift_ref = spectrum.read_hdf(infile['shift_ref'])

        if 'shift_data' in infile:
            sm.shift_data = h5plus.read_dict(infile['shift_data'])

        # ------------------------ Match results ------------------------
        if 'regions' in infile:
            regions = infile['regions'][:]
            sm.regions = []
            for i in range(regions.shape[0]):
                sm.regions.append((regions[i, 0], regions[i, 1]))

        if 'match_results' in infile:
            match_results = pd.DataFrame.from_records(
                infile['match_results'][:], index='lib_index')
            # decode strings
            for (col_name, dt) in match_results.dtypes.iteritems():
                if dt == 'object':
                    match_results[col_name] = match_results[col_name]\
                        .str.decode('utf-8')
            match_results['lib_index'] = match_results.index
            sm.match_results = match_results

        # ------------------------ Match results ------------------------
        if 'num_best' in infile:
            sm.num_best = infile['num_best'].value

        if 'coeffs' in infile:
            sm.coeffs = infile['coeffs'][:]

        if 'lincomb_regions' in infile:
            lincomb_regions = infile['lincomb_regions'][:]
            sm.lincomb_regions = []
            for i in range(lincomb_regions.shape[0]):
                sm.lincomb_regions.append((lincomb_regions[i, 0],
                                           lincomb_regions[i, 1]))

        if 'lincomb' in infile and sm.lincomb_regions is not None:
            grp = infile['lincomb']
            ref_idxs = []
            lincomb_matches = []
            lincomb_results = []
            for i in range(len(sm.lincomb_regions)):
                reg = sm.lincomb_regions[i]
                subgrp = grp['{0:.0f}'.format(reg[0])]
                ref_idxs.append(subgrp['ref_idxs'][:])

                mt = match.MatchLincomb.read_hdf(subgrp)
                lincomb_matches.append(mt)

                res_grp = subgrp['results']
                reg_results = {}
                for k in res_grp:
                    reg_results[k] = res_grp[k].value
                lincomb_results.append(reg_results)

            sm.ref_idxs = ref_idxs
            sm.lincomb_matches = lincomb_matches
            sm.lincomb_results = lincomb_results

        if 'results_nodetrend' in infile:
            res_grp = infile['results_nodetrend']
            results_nodetrend = {}
            for k in res_grp:
                results_nodetrend[k] = res_grp[k].value

            sm.results_nodetrend = results_nodetrend

        if 'results' in infile:
            res_grp = infile['results']
            results = {}
            for k in res_grp:
                results[k] = res_grp[k].value

            sm.results = results

        if is_path:
            infile.close()

        return sm

    def to_hdf(self, outfile):
        """Saves the current state of the SpecMatch object to an hdf file.

        Args:
            outfile (str or h5 file): Output path or file handle.
        """
        print("Saving SpecMatch object to HDF file")
        # Allow either a string or h5 file object ot be passed.
        is_path = False
        if isinstance(outfile, string_types):
            outfile = h5py.File(outfile, 'w')
            is_path = True

        # Save target spectrum
        if self.target is not None:
            outfile.create_group('target')
            self.target.to_hdf(outfile['target'])

        # ------------------------ Shift results ------------------------
        # Save unshifted spectrum
        if self._shifted is True:
            # (Check if unshifted and shifted targets are different)
            outfile.create_group('unshifted')
            self.target_unshifted.to_hdf(outfile['unshifted'])

        # Save shift reference
        if self.shift_ref is not None:
            outfile.create_group('shift_ref')
            self.shift_ref.to_hdf(outfile['shift_ref'])

        # Save shift data
        if len(self.shift_data) > 0:
            grp = outfile.create_group('shift_data')
            for k, v in self.shift_data.items():
                grp[k] = v

        # ------------------------ Match results ------------------------
        # Wavelength regions
        if self.regions is not None:
            outfile['regions'] = self.regions

        # Save match results by converting to record array
        if not self.match_results.empty:
            if 'lib_index' in self.match_results:
                self.match_results.drop('lib_index', inplace=True, axis=1)
            match_rec = self.match_results.to_records()
            dt = match_rec.dtype.descr

            for i in range(len(dt)):
                if dt[i][1] == "|O":
                    # max string length = 5000
                    dt[i] = (str(dt[i][0]), 'S5000')
                else:
                    dt[i] = (str(dt[i][0]), dt[i][1])

            match_rec = np.array(match_rec, dtype=dt)

            outfile['match_results'] = match_rec

            # restore lib_index column
            self.match_results['lib_index'] = self.match_results.index

        # ----------------------- Lincomb results -----------------------
        if self.num_best is not None:
            outfile['num_best'] = self.num_best
        if self.coeffs is not None:
            outfile['coeffs'] = self.coeffs
        if self.lincomb_regions is not None:
            outfile['lincomb_regions'] = self.lincomb_regions
        if len(self.lincomb_matches) > 0:
            grp = outfile.create_group('lincomb')
            for i in range(len(self.lincomb_regions)):
                reg = self.lincomb_regions[i]
                subgrp = grp.create_group('{0:.0f}'.format(reg[0]))
                subgrp['ref_idxs'] = self.ref_idxs[i]
                self.lincomb_matches[i].to_hdf(subgrp)
                res_grp = subgrp.create_group('results')
                for k, v in self.lincomb_results[i].items():
                    res_grp[k] = v

        # Save averaged results
        if len(self.results_nodetrend) > 0:
            res_grp = outfile.create_group('results_nodetrend')
            for k, v in self.results_nodetrend.items():
                res_grp[k] = v

        if len(self.results) > 0:
            res_grp = outfile.create_group('results')
            for k, v in self.results.items():
                res_grp[k] = v

        if is_path:
            outfile.close()

    def results_to_txt(self, outfile, verbose=False):
        """Save final results to text.

        Args:
            outfile (str or file object): Output file
            verbose (bool): Whether to only print the final derived results
                or to also include not-detrended results and list of best
                matching spectra for each region
        """
        if isinstance(outfile, string_types):
            f = open(outfile, 'w')
        else:
            f = outfile

        f.write('Derived Parameters\n')
        f.write('------------------\n')
        f.write('Teff: {0:.0f} +/- {1:.0f} K\n'.format(
            self.results['Teff'], self.results['u_Teff']))
        f.write('Radius: {0:.3f} +/- {1:.3f} Rsun\n'.format(
            self.results['radius'], self.results['u_radius']))
        f.write('[Fe/H]: {0:.2f} +/- {1:.2f} dex\n'.format(
            self.results['feh'], self.results['u_feh']))

        if verbose:
            f.write('\n')
            f.write('Parameters before detrending\n')
            f.write('----------------------------\n')
            f.write('Teff: {0:.0f} +/- {1:.0f} K\n'.format(
                self.results_nodetrend['Teff'],
                self.results_nodetrend['u_Teff']))
            f.write('Radius: {0:.3f} +/- {1:.3f} Rsun\n'.format(
                self.results_nodetrend['radius'],
                self.results_nodetrend['u_radius']))
            f.write('[Fe/H]: {0:.2f} +/- {1:.2f} dex\n'.format(
                self.results_nodetrend['feh'],
                self.results_nodetrend['u_feh']))

            f.write('\n')
            f.write('Best Matching Spectra\n')
            f.write('---------------------\n')
            for i in range(len(self.regions)):
                f.write('Region {0}:\n'.format(self.regions[i]))
                mt = self.lincomb_matches[i]
                for j in range(mt.num_refs):
                    ref = mt.refs[j]
                    f.write('\t#{0:d}: {1}, '.format(j, ref.name))
                    f.write('chi^2 = {0:.3f}, '.format(mt.ref_chisq[j]))
                    f.write('c_{0:d} = {1:.3f}\n'.format(j, mt.coeffs[j]))
                f.write('Final chi^2 = {0:.3f}\n'.format(mt.best_chisq))

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

        # ------------------------ Shift results ------------------------
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
        if len(self.shift_data) > 0:
            shift_data = self.shift_data.copy()
            col_list = []
            num_orders = shift_data.pop('num_orders')
            col_list.append(fits.Column(name='num_orders', format='J',
                                        array=[num_orders]))

            num_sects = []
            for i in range(num_orders):
                num_sects.append(shift_data.pop('order_{0:d}/num_sections'
                                                .format(i)))
            col_list.append(fits.Column(name='num_sects', format='J',
                                        array=num_sects))
            n_sec = max(num_sects)

            for k in ['center_pix', 'lag', 'fit']:
                col_list.append(fits.Column(name=k,
                                format='{0:d}E'.format(n_sec),
                                array=shift_data.pop(k)))

            # Save individual fit data
            for k in shift_data.keys():
                col_list.append(fits.Column(name=k, format='D',
                                            array=shift_data[k]))

            shift_hdu = fits.BinTableHDU.from_columns(col_list)
            shift_hdu.name = 'SHIFTDATA'
            hdulist.append(shift_hdu)

        # ------------------------ Match results ------------------------
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

        # ----------------------- Lincomb results -----------------------

    # ----------------------------- Shift plots ----------------------------- #
    def plot_shifted_spectrum(self, wavlim=(5158, 5172)):
        """Plot a comparison of the shifted, reference and unshifted spectra.

        Args:
            wavlim (tuple or list of tuples, optional):
                Wavelength limits to plot.
        """
        if not isinstance(wavlim, list) and not isinstance(wavlim, np.ndarray):
            wavlim = [wavlim]

        num_regions = np.shape(wavlim)[0]

        for i in range(num_regions):
            plt.subplot(num_regions, 1, i + 1)

            target = self.target.cut(*wavlim[i])

            targid = target.name
            target.plot(offset=1, text='Target (shifted): {0}'.format(targid))

            target_unshifted = self.target_unshifted.cut(*wavlim[i])
            target_unshifted.plot(offset=0, normalize=True,
                                  text='Target (unshifted)',
                                  plt_kw={'color': 'forestgreen'})

            if self.shift_ref is not None:
                # Plot shift reference
                shift_ref = self.shift_ref.cut(*wavlim[i])
                refid = shift_ref.name
                shift_ref.plot(offset=1.5, text='Reference: {0}'.format(refid),
                               plt_kw={'color': 'firebrick'})

                # Plot residuals
                if (target.w[0] > shift_ref.w[0]) \
                        or (target.w[-1] < shift_ref.w[-1]):
                    target = target.extend(shift_ref.w)

                if len(shift_ref.w) < len(target.w):
                    target = target.cut(shift_ref.w[0], shift_ref.w[-1])
                plt.plot(target.w, shift_ref.s - target.s, '-', color='purple')
                plots.annotate_spectrum('Residuals', spec_offset=-1)

            if self.lib.nso is not None and self.shift_ref.name != 'NSO':
                nso = self.lib.nso.cut(*wavlim[i])
                nso.plot(offset=2, text='NSO', plt_kw={'color': 'c'})
                plt.ylim(-0.5, 3.5)
            else:
                plt.ylim(-0.5, 2.7)

            plt.xlim(wavlim[i])

    def plot_shift_lags(self, orders='all', legend=True):
        """Plot the lags for each order as a function of pixel number.

        Args:
            orders (str or int or list of ints):
                'all' (default): Plot all orders
                int: Plot only given order
                list of ints: Plot the orders provided.
        """
        lags = self.shift_data['lag']
        center_pix = self.shift_data['center_pix']
        fit = self.shift_data['fit']

        total_orders = lags.shape[0]
        if isinstance(orders, int):
            orders = [orders]
        if orders == 'all':
            orders = np.arange(total_orders)

        colormap = plt.cm.nipy_spectral
        num_orders = len(orders)
        for i in range(num_orders):
            color = colormap(0.9 * i / num_orders + 0.1)
            order = orders[i]
            plt.plot(
                center_pix[order], lags[order], '_', ms=6, mew=2, color=color
            )
            plt.plot(
                center_pix[order], fit[order], '-', color=color,
                label='{0:d}'.format(order)
            )

        plt.xlabel('Pixel number')
        plt.ylabel('Shift (pixels)')
        if legend:
            plt.legend(loc='best', ncol=2, fontsize='small',title='Order')

    def plot_xcorr(self, order=0, highlightpeak=False, legend=True):
        """Plot the correlation array for each section of the given order.

        Args:
            order (int): Order on chip to plot. Defaults to 0.
            highlightpeak (bool): Whether to highlight the peak value.
        """
        num_sects = self.shift_data['order_{0:d}/num_sections'.format(order)]
        colormap = plt.cm.nipy_spectral

        for i in range(num_sects):
            color = colormap(0.9 * i / num_sects + 0.1)

            xcorr = self.shift_data['order_{0:d}/sect_{1:d}/xcorr'
                                    .format(order, i)]
            lag_arr = self.shift_data['order_{0:d}/sect_{1:d}/lag_arr'
                                      .format(order, i)]
            plt.plot(lag_arr, xcorr, label="{0:d}".format(i),color=color)

            if highlightpeak and len(xcorr) > 0:
                max_corr = np.argmax(xcorr)
                mew = plt.rcParams['lines.markeredgewidth'] * 2
                plt.plot(
                    lag_arr[max_corr], xcorr[max_corr],'x',color=color, 
                    mew=mew,zorder=10
                )

        if legend:
            plt.legend(loc='upper left', fontsize='small',title='Section')

        plt.xlim(-600, 600)
        plt.xlabel('Shift (pixels)')
        plt.ylabel('Correlation')

    # ----------------------------- Match plots ----------------------------- #
    def plot_chi_squared_surface(self, region=0, num_best=None):
        """Plot the chi-squared surface from the pairwise matching procedure.

        Creates a three-column plot of the best chi-squared obtained with
        each library spectrum as a function of the library parameters.

        Args:
            region (int or tuple): The wavelength region to plot.
                If an int is given, refers to the index of the region in
                `self.regions`.
                If a tuple is given, should be present in `self.regions`.
            num_best (optional [int]): Number of best spectra to highlight.
                (default: `self.num_best`)
        """
        if self.match_results.empty:
            return

        if isinstance(region, int):
            region = self.regions[region]
        elif isinstance(region, tuple) and region not in self.regions:
            raise ValueError("Region {0} was not a region.".format(region))

        if len(self.regions) == 1:
            cs_col = 'chi_squared'
        else:
            cs_col = 'chi_squared_{0:.0f}'.format(region[0])
        self.match_results.sort_values(by=cs_col, inplace=True)

        if num_best is None:
            num_best = self.num_best

        ax1 = plt.subplot(131)
        plt.semilogy()
        plt.plot(self.match_results['Teff'], self.match_results[cs_col], '.',
                 color='blue')
        plt.plot(self.match_results['Teff'].head(num_best),
                 self.match_results[cs_col].head(num_best), 'r.')
        plt.ylabel(r'$\chi^2$')
        plt.xlabel(r'$T_{eff}$ (K)')
        plots.label_axes(param_x='Teff')

        plt.subplot(132, sharey=ax1)
        plt.semilogy()
        plt.plot(self.match_results['radius'], self.match_results[cs_col], '.',
                 color='blue')
        plt.plot(self.match_results['radius'].head(num_best),
                 self.match_results[cs_col].head(num_best), 'r.')
        plots.label_axes(param_x='radius')

        plt.subplot(133, sharey=ax1)
        plt.semilogy()
        plt.plot(self.match_results['feh'], self.match_results[cs_col], '.',
                 color='blue')
        plt.plot(self.match_results['feh'].head(num_best),
                 self.match_results[cs_col].head(num_best), 'r.')
        plots.label_axes(param_x='feh')

    def plot_best_match_spectra(self, region=0, wavlim='all', num_best=None):
        """Plots the reference, modified reference and residuals for each of the
        best matches.

        Args:
            region (int or tuple): The region used in matching.
                If an int is given, refers to the index of the region in
                `self.regions`.
                If a tuple is given, should be present in `self.regions`.
            wavlim (str or tuple): The wavelength limits to plot.
                If 'all' is given, plots the entire `region`.
            num_best (optional [int]): Number of best spectra to highlight.
                (default: `self.num_best`)
        """
        if self.match_results.empty:
            return

        # get region
        if isinstance(region, int):
            region = self.regions[region]
        elif isinstance(region, tuple) and region not in self.regions:
            raise ValueError("Region {0} was not a region.".format(region))

        if len(self.regions) == 1:
            fit_col = 'fit_params'
            cs_col = 'chi_squared'
        else:
            fit_col = 'fit_params_{0:.0f}'.format(region[0])
            cs_col = 'chi_squared_{0:.0f}'.format(region[0])

        self.match_results.sort_values(by=cs_col, inplace=True)

        # get wavlength region to plot
        if isinstance(wavlim, tuple):
            # cannot exceed region
            plotwavlim = (max(wavlim[0], region[0]), min(wavlim[1], region[1]))
            if plotwavlim[0] >= plotwavlim[1]:
                raise ValueError("Wavlim {0} is not within region {1}".format(
                    wavlim, region))
        else:
            plotwavlim = region

        # get num_best
        if num_best is None:
            num_best = self.num_best

        shared_ax = None
        for i in range(num_best):
            if shared_ax is None:
                ax = plt.subplot(num_best, 1, i + 1)
                shared_ax = ax
            else:
                ax = plt.subplot(num_best, 1, i + 1, sharex=shared_ax)
            ref_idx = self.match_results.iloc[i]['lib_index']
            ref = self.lib.get_spectrum(ref_idx)
            params = lmfit.Parameters()
            params.loads(self.match_results.iloc[i][fit_col])
            mt = match.Match(self.target.cut(*region), ref.cut(*region))
            mt.load_params(params)
            mt.plot()
            if i != num_best - 1:
                plt.setp(ax.get_xticklabels(), visible=False)
            ax.set_xlabel('')
            ax.set_ylabel('')
            plt.xlim(plotwavlim)

        # Label
        fig = plt.gcf()
        axes = fig.axes
        axes[-1].set_xlabel('Wavelength (Angstroms)', fontsize='large')
        fig.text(0.05, 0.5, 'Normalized Flux (Arbitrary Offset)',
                 ha='center', va='center', rotation='vertical',
                 fontsize='large')
        plt.tight_layout(rect=[0.05, 0, 1, 0.95])

    def plot_references(self, region=0, num_best=None, verbose=True):
        """Plots the location of the best references used in the linear
        combination step.

        Args:
            region (int or tuple): The region used in matching.
                If an int is given, refers to the index of the region in
                `self.regions`.
                If a tuple is given, should be present in `self.regions`.
            num_best (optional [int]): Number of best spectra to highlight.
                (default: `self.num_best`)
            verbose (optional [bool]): Whether to annotate the points with
                the lincomb coefficients
        """
        if self.match_results is None:
            return

        # get region
        if isinstance(region, int):
            region_num = region
            region = self.regions[region_num]
        elif isinstance(region, tuple):
            if region in self.regions:
                region_num = self.regions.index(region)
            else:
                # Raise exception of region was not found.
                raise ValueError("Region {0} was not a region.".format(region))

        if len(self.regions) == 1:
            cs_col = 'chi_squared'
        else:
            cs_col = 'chi_squared_{0:.0f}'.format(region[0])
        self.match_results.sort_values(by=cs_col, inplace=True)

        if num_best is None:
            num_best = self.num_best

        def _plot_ref_params(paramx, paramy, zoom=False):
            plt.plot(self.match_results[paramx], self.match_results[paramy],
                     'b.', alpha=0.6, label='_nolegend_')
            plt.plot(self.match_results.head(num_best)[paramx],
                     self.match_results.head(num_best)[paramy], '^',
                     color='forestgreen', label='Best Matches')
            if zoom:
                plots.set_tight_lims(self.match_results.head(num_best)[paramx],
                                     self.match_results.head(num_best)[paramy])

                if verbose and self.coeffs is not None:
                    for i in range(self.num_best):
                        p = self.match_results.iloc[i]
                        plots.annotate_point(p[paramx], p[paramy],
                                             '{0:.3f}'.format(
                                             self.coeffs[region_num][i]))

        gs = gridspec.GridSpec(2, 2)
        plt.subplot(gs[0])
        _plot_ref_params('Teff', 'radius')
        if len(self.results) > 0:
            plt.plot(self.lincomb_results[region_num]['Teff'],
                     self.lincomb_results[region_num]['radius'], 's',
                     color='purple', label='Derived Parameters')
        plots.label_axes('Teff', 'radius')

        plt.subplot(gs[1])
        _plot_ref_params('Teff', 'radius', zoom=True)
        if len(self.results) > 0:
            plt.plot(self.lincomb_results[region_num]['Teff'],
                     self.lincomb_results[region_num]['radius'], 's',
                     color='purple', label='Derived Parameters')
        plots.label_axes('Teff', 'radius', rescale=False)

        plt.subplot(gs[2])
        _plot_ref_params('feh', 'radius')
        if len(self.results) > 0:
            plt.plot(self.lincomb_results[region_num]['feh'],
                     self.lincomb_results[region_num]['radius'], 's',
                     color='purple', label='Derived Parameters')
        plots.label_axes('feh', 'radius')

        plt.subplot(gs[3])
        _plot_ref_params('feh', 'radius', zoom=True)
        if len(self.results) > 0:
            plt.plot(self.lincomb_results[region_num]['feh'],
                     self.lincomb_results[region_num]['radius'], 's',
                     color='purple', label='Derived Parameters')
        plots.label_axes('feh', 'radius', rescale=False)

    def plot_lincomb(self, region=0, wavlim='all'):
        """Plots the spectra used to synthesize the linear combination,
        the synthesized spectrum, and residuals.
        """
        if len(self.lincomb_matches) == 0:
            return

        # get region
        if isinstance(region, int):
            region_num = region
            region = self.lincomb_regions[region_num]
        elif isinstance(region, tuple):
            if region in self.lincomb_regions:
                region_num = self.lincomb_regions.index(region)
            else:
                # Raise exception of region was not found.
                raise ValueError("Region {0} was not a region.".format(region))

        # get wavlength region to plot
        if isinstance(wavlim, tuple):
            # cannot exceed region
            plotwavlim = (max(wavlim[0], region[0]), min(wavlim[1], region[1]))
            if plotwavlim[0] >= plotwavlim[1]:
                raise ValueError("Wavlim {0} is not within region {1}".format(
                    wavlim, region))
        else:
            plotwavlim = region

        mt_lincomb = self.lincomb_matches[region_num]

        mt_lincomb.plot()
        plt.xlim(plotwavlim)
