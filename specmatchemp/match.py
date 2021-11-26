"""
@filename match.py

Defines the Match class
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import lmfit
from scipy.interpolate import LSQUnivariateSpline
from scipy.ndimage.filters import convolve1d

import specmatchemp.kernels
from specmatchemp import spectrum
from specmatchemp.spectrum import Spectrum
from specmatchemp import plots


class Match(object):
    """The Match class used for matching two spectra

    Attributes:
        target (Spectrum): Target spectrum
        targ_mod (Spectrum): Modified target spectrum, broadened to account
            for negative relative vsini
        reference (Spectrum): Reference spectrum
        ref_mod (Spectrum): Modified reference, after broadening and fitting a
            spline to correct for continuum differences
        spl (np.ndarray): Spline used to fit continuum level
        best_params (lmfit.Parameters): Parameters used to create the best
            match.
        best_chisq (float): Chi-squared from the best match

    Args:
        target (Spectrum): Target spectrum
        reference (Spectrum): Reference spectrum
        mode: default (unnormalized chi-square),
              normalized (normalized chi-square)
        opt: lm (Levenberg-Marquadt optimization), nelder (Nelder-Mead)
    """
    def __init__(self, target, reference, mode='default', opt='nelder'):
        # Ensure target, reference are on common wavelength scale
        if not np.allclose(target.w, reference.w):
            print("Target and reference are on different wavelength scales.")
            raise ValueError
        # common wavelength scale
        self.w = np.copy(target.w)

        # target, reference and modified spectra
        self.target = target.copy()
        self.targ_mod = target.copy()
        self.reference = reference.copy()
        self.ref_mod = reference.copy()

        # replace nans with continuum
        self.target.s[np.isnan(self.target.s)] = 1.
        self.target.serr[np.isnan(self.target.serr)] = 1.
        self.reference.s[np.isnan(self.reference.s)] = 1.
        self.reference.serr[np.isnan(self.reference.serr)] = 1.

        self.best_params = lmfit.Parameters()
        self.best_chisq = np.NaN
        self.mode = mode
        self.opt = opt

        # add spline knots
        num_knots = 5
        interval = int(len(self.w)/(num_knots+1))
        # Add spline positions
        self.knot_x = []
        for i in range(1, num_knots+1):
            self.knot_x.append(self.w[interval*i])
        self.knot_x = np.array(self.knot_x)

    def create_model(self, params):
        """
        Creates a tweaked model based on the parameters passed,
        based on the reference spectrum.
        Stores the tweaked model in spectra.s_mod and serr_mod.
        """
        self.targ_mod.s = np.copy(self.target.s)
        self.targ_mod.serr = np.copy(self.target.serr)
        self.ref_mod.s = np.copy(self.reference.s)
        self.ref_mod.serr = np.copy(self.reference.serr)

        # Apply broadening kernel
        vsini = params['vsini'].value

        if vsini >= 1.0:
            self.ref_mod = self.broaden(vsini, self.ref_mod)
        elif vsini < 1.0:
            # Threshold at 1 because kernel does nothing for abs(vsini) < 1.0
            self.targ_mod = self.broaden(2.0 - vsini, self.targ_mod)

        # Use linear least squares to fit a spline
        spline = LSQUnivariateSpline(self.w, self.targ_mod.s / self.ref_mod.s,
                                     self.knot_x)
        self.spl = spline(self.w)

        self.ref_mod.s *= self.spl
        self.ref_mod.serr *= self.spl

    def load_params(self, params):
        """
        Method to create a model based on pre-determined parameters,
        storing it as the best fit model.
        """
        self.best_chisq = self.objective(params)
        self.best_params = params

    def broaden(self, vsini, spec):
        """
        Applies a broadening kernel to the given spectrum (or error)

        Args:
            vsini (float): vsini to determine width of broadening
            spec (Spectrum): spectrum to broaden
        Returns:
            broadened (Spectrum): Broadened spectrum
        """
        SPEED_OF_LIGHT = 2.99792e5
        dv = (self.w[1]-self.w[0])/self.w[0]*SPEED_OF_LIGHT
        n = 151     # fixed number of points in the kernel
        varr, kernel = specmatchemp.kernels.rot(n, dv, vsini)
        # broadened = signal.fftconvolve(spec, kernel, mode='same')

        spec.s = convolve1d(spec.s, kernel)
        spec.serr = convolve1d(spec.serr, kernel)

        return spec

    def objective(self, params):
        """
        Objective function evaluating goodness of fit given the passed
        parameters.

        Args:
            params
        Returns:
            Reduced chi-squared value between the target spectra and the
            model spectrum generated by the parameters
        """
        self.create_model(params)

        # Calculate residuals (data - model)
        if self.mode == 'normalized':
            residuals = ((self.targ_mod.s - self.ref_mod.s) /
                         np.sqrt(self.targ_mod.serr**2 + self.ref_mod.serr**2))
        else:
            residuals = (self.targ_mod.s - self.ref_mod.s)

        chi_square = np.sum(residuals**2)

        if self.opt == 'lm':
            return residuals
        elif self.opt == 'nelder':
            return chi_square

    def best_fit(self, params=None, allow_negvsini=False):
        """
        Calculates the best fit model by minimizing over the parameters:
        - spline fitting to the continuum
        - rotational broadening
        """
        if params is None:
            params = lmfit.Parameters()

        # Rotational broadening parameters
        if allow_negvsini:
            params.add('vsini', value=1.0, min=-5.0, max=10.0)
        else:
            params.add('vsini', value=1.0, min=1.0, max=10.0)

        # Spline parameters
        params = add_spline_positions(params, self.knot_x)

        # Perform fit
        if self.opt == 'lm':
            out = lmfit.minimize(self.objective, params)
            self.best_chisq = np.sum(self.objective(out.params)**2)
        elif self.opt == 'nelder':
            out = lmfit.minimize(self.objective, params, method='nelder')
            self.best_chisq = self.objective(out.params)

        self.best_params = out.params

        return self.best_chisq

    def best_residuals(self):
        """Returns the residuals between the target spectrum and best-fit
        spectrum.

        Returns:
            np.ndarray
        """
        if self.mode == 'normalized':
            return ((self.targ_mod.s - self.ref_mod.s) /
                    np.sqrt(self.targ_mod.serr**2 + self.ref_mod.serr**2))
        else:
            return (self.targ_mod.s - self.ref_mod.s)  # data - model

    def get_spline_positions(self):
        """Wrapper function for getting spline positions

        Returns:
            knotx (np.ndarray)
        """
        return get_spline_positions(self.best_params)

    def plot(self, offset=1, verbose=True):
        if verbose:
            vsini = get_physical_vsini(self.best_params['vsini'].value)
            labels = {'target': 'Target: {0}'.format(self.targ_mod.name),
                    'reference': 'Reference: {0}'.format(self.reference.name),
                    'modified': r'Reference (modified): $v\sin i = {0:.2f}$'
                                .format(vsini),
                    'residuals': r'Residuals: $\chi^2 = {0:.3f}$'
                                .format(self.best_chisq)}
        else:
            labels = {'target': 'Target',
                      'reference': 'Reference',
                      'modified': 'Reference (Modified)',
                      'residuals': 'Residuals'}

        self.targ_mod.plot(text=labels['target'],
                           plt_kw={'color': 'royalblue'})
        self.ref_mod.plot(offset=1*offset, plt_kw={'color': 'forestgreen'},
                          text=labels['modified'])
        self.reference.plot(offset=2*offset, plt_kw={'color': 'firebrick'},
                            text=labels['reference'])

        plt.plot(self.targ_mod.w, self.ref_mod.s-self.targ_mod.s,
                 '-', color='black')
        plots.annotate_spectrum(labels['residuals'], spec_offset=-1)

    def to_hdf(self, outfile):
        """Saves the Match object to an hdf file.

        Args:
            outfile (str or h5 file): Output path or file handle.
        """
        # Allow either a string or h5 file object to be passed.
        is_path = False
        if isinstance(outfile, str):
            outfile = h5py.File(outfile, 'w')
            is_path = True

        # Specmatch-Emp version
        outfile['smemp-version'] = specmatchemp.SPECMATCH_VERSION.lstrip('v')

        # Save target
        outfile.create_group('target')
        self.target.to_hdf(outfile['target'])

        # Save modified target
        outfile.create_group('targ_mod')
        self.targ_mod.to_hdf(outfile['targ_mod'])

        # Save reference
        outfile.create_group('reference')
        self.reference.to_hdf(outfile['reference'])

        # Save modified
        outfile.create_group('ref_mod')
        self.ref_mod.to_hdf(outfile['ref_mod'])

        # Save best-fit parameters
        outfile['best_params'] = self.best_params.dumps()
        outfile['spline'] = self.spl
        outfile['best_chisq'] = self.best_chisq

        if is_path:
            outfile.close()

    @classmethod
    def read_hdf(cls, infile):
        """Reads the Match object from an HDF file.

        Args:
            infile (str or h5 file): Input path or file handle.
        """
        # Allow either a string or h5 file object to be passed.
        is_path = False
        if isinstance(infile, str):
            infile = h5py.File(infile, 'r')
            is_path = True

        # Backward compatibility
        if 'smemp-version' not in infile or infile['smemp-version'] < 0.4:
            # Read target
            target = spectrum.read_hdf(infile['target'])
            targ_mod = target.copy()
            # Read reference
            reference = spectrum.read_hdf(infile['reference'])

            # Read modified
            ref_mod = spectrum.read_hdf(infile['modified'])
        else:
            target = spectrum.read_hdf(infile['target'])
            targ_mod = spectrum.read_hdf(infile['targ_mod'])
            reference = spectrum.read_hdf(infile['reference'])
            ref_mod = spectrum.read_hdf(infile['ref_mod'])

        # Read best-fit parameters
        best_params = lmfit.Parameters()
        best_params.loads(infile['best_params'].value)
        spl = infile['spline'][:]
        best_chisq = infile['best_chisq'].value

        mt = cls(target, reference)
        mt.load_params(best_params)

        if not np.allclose(ref_mod.s, mt.ref_mod.s) \
                or not np.allclose(targ_mod.s, mt.targ_mod.s) \
                or not np.allclose(spl, mt.spl) \
                or mt.best_chisq != best_chisq:
            print("Warning: Saved model and recreated model are not " +
                  "identical. Model may have been created by a different " +
                  "version of SpecMatch-Emp.")

        if is_path:
            infile.close()

        return mt


class MatchLincomb(Match):
    """MatchLincomb class used to match a linear combination of spectra to
    match a target spectrum.

    Attributes:
        target (Spectrum): Target spectrum
        refs (list of Spectrum): List of reference spectra
        vsini (np.ndarray): Array of vsini used to broaden the reference
            spectra.
        spl (np.ndarray): Spline used for continuum fitting.
        modified (Spectrum): Best linear combination of spectra.
        best_params (lmfit.Parameters): Parameters used to create the best
            match.
        best_chisq (float): Chi-squared from the best match.

    Args:
        target (Spectrum): Target spectrum
        refs (list of Spectrum): Array of reference spectra
        vsini (np.ndarray): array containing vsini broadening for each
                            reference spectrum
    """
    def __init__(self, target, refs, vsini, mode='default', ref_chisq=None):
        # Ensure all references and target are on the same wavelength scale
        for i in range(len(refs)):
            if not np.allclose(target.w, refs[i].w):
                print("Target and reference {0:d} are on different".format(i) +
                      "wavelength scales.")
                raise ValueError

        self.w = np.copy(target.w)
        self.target = target.copy()
        self.targ_mod = target.copy()

        # replace nans with continuum
        self.targ_mod.s[np.isnan(self.targ_mod.s)] = 1.
        self.targ_mod.serr[np.isnan(self.targ_mod.serr)] = 1.

        self.num_refs = len(refs)
        self.refs = []
        for i in range(self.num_refs):
            self.refs.append(refs[i].copy())
            self.refs[i].s[np.isnan(self.refs[i].s)] = 1.
            self.refs[i].serr[np.isnan(self.refs[i].serr)] = 1.

        self.ref_chisq = ref_chisq

        self.original_vsini = vsini.copy()
        self.vsini = vsini
        self.phys_vsini = [get_physical_vsini(v) for v in vsini]

        # Broaden target spectrum to make all vsini positive
        min_vsini = np.min(self.phys_vsini)
        if min_vsini < 0.0:
            self.targ_mod = self.broaden(-min_vsini, self.targ_mod)
            for idx, v in enumerate(self.vsini):
                if v >= 1.0:
                    self.vsini[idx] = v - min_vsini
                else:
                    # Correction for discontinuity
                    self.vsini[idx] = v - min_vsini - 1.0

        self.refs_broadened = []
        for i in range(self.num_refs):
            self.refs_broadened.append(self.refs[i].copy())
            if vsini[i] >= 1.0:
                self.refs_broadened[i] = self.broaden(vsini[i],
                                                      self.refs_broadened[i])

        self.ref_mod = Spectrum(self.w, np.zeros_like(self.w),
                                name='Linear Combination {0:d}'
                                .format(self.num_refs))

        self.best_params = lmfit.Parameters()
        self.best_chisq = np.NaN
        self.mode = mode
        self.opt = 'nelder'

        # add spline knots
        num_knots = 5
        interval = int(len(self.w)/(num_knots+1))
        # Add spline positions
        self.knot_x = []
        for i in range(1, num_knots+1):
            self.knot_x.append(self.w[interval*i])
        self.knot_x = np.array(self.knot_x)

    def create_model(self, params):
        """
        Creates a tweaked model based on the parameters passed,
        based on the reference spectrum.
        Stores the tweaked model in spectra.s_mod and serr_mod.
        """
        self.ref_mod.s = np.zeros_like(self.w)
        self.ref_mod.serr = np.zeros_like(self.w)

        # create the model from a linear combination of the reference spectra
        coeffs = get_lincomb_coeffs(params)

        for i in range(self.num_refs):
            self.ref_mod.s += self.refs_broadened[i].s * coeffs[i]
            self.ref_mod.serr += self.refs_broadened[i].serr * coeffs[i]

        # Use linear least squares to fit a spline
        spline = LSQUnivariateSpline(self.w, self.targ_mod.s / self.ref_mod.s,
                                     self.knot_x)
        self.spl = spline(self.w)

        self.ref_mod.s *= self.spl
        self.ref_mod.serr *= self.spl

    def objective(self, params):
        """Objective function evaluating goodness of fit given the passed
        parameters.

        Args:
            params
        Returns:
            Reduced chi-squared value between the target spectra and the
            model spectrum generated by the parameters
        """
        chi_square = super(MatchLincomb, self).objective(params)

        # Add a Gaussian prior
        sum_coeff = np.sum(get_lincomb_coeffs(params))

        WIDTH = 1e-2
        chi_square += (sum_coeff - 1)**2 / (2 * WIDTH**2)

        return chi_square

    def best_fit(self):
        """
        Calculates the best fit model by minimizing over the parameters:
        - Coefficients of reference spectra
        - spline fitting to the continuum
        - rotational broadening
        """
        params = lmfit.Parameters()

        # Linear combination parameters
        params = add_lincomb_coeffs(params, self.num_refs)

        # Spline parameters
        params = add_spline_positions(params, self.knot_x)

        # vsini
        params = add_vsini(params, self.vsini)

        # Minimize chi-squared
        out = lmfit.minimize(self.objective, params, method='nelder')

        # Save best fit parameters
        self.best_params = out.params
        self.best_chisq = self.objective(self.best_params)
        self.coeffs = self.get_lincomb_coeffs()

        return self.best_chisq

    def get_vsini(self):
        """Wrapper function to get vsini list from MatchLincomb object

        Returns:
            vsini (np.ndarray)
        """
        return get_vsini(self.best_params)

    def get_lincomb_coeffs(self):
        """Wrapper function to get lincomb coefficients from MatchLincomb object

        Returns:
            coeffs (np.ndarray)
        """
        return get_lincomb_coeffs(self.best_params)

    def plot(self, verbose=True):
        # create labels
        if verbose:
            labels = {'target': 'Target: {0}'.format(self.target.name),
                      'modified': r'Linear Combination',
                      'residuals': r'Residuals: $\chi^2 = {0:.3f}$'
                                   .format(self.best_chisq)}

            coeffs = self.get_lincomb_coeffs()

            for i in range(self.num_refs):
                if self.ref_chisq is None:
                    labels['ref_{0:d}'.format(i)] = (
                        'Reference: {0}, '.format(self.refs[i].name) +
                        r'$v\sin i = {0:.2f}$, '.format(self.phys_vsini[i]) +
                        r'$c_{0:d} = {1:.3f}$'.format(i, coeffs[i]))

                else:
                    labels['ref_{0:d}'.format(i)] = (
                        'Reference: {0}, '.format(self.refs[i].name) +
                        r'$v\sin i = {0:.2f}$, '.format(self.phys_vsini[i]) +
                        r'$\chi^2 = {0:.2f}$, '.format(self.ref_chisq[i]) +
                        r'$c_{0:d} = {1:.3f}$'.format(i, coeffs[i]))
        else:
            labels = {'target': 'Target', 'modified': 'Reference (Modified)',
                      'residuals': 'Residuals'}

            for i in range(self.num_refs):
                labels['ref_{0:d}'.format(i)] = 'Reference {0:d}'.format(i)

        # Plot spectra
        self.targ_mod.plot(plt_kw={'color': 'royalblue'},
                           text=labels['target'])
        self.ref_mod.plot(offset=0.5, plt_kw={'color': 'forestgreen'},
                          text=labels['modified'])

        for i in range(self.num_refs):
            self.refs[i].plot(offset=1.5+i*0.5, plt_kw={'color': 'firebrick'},
                              text=labels['ref_{0:d}'.format(i)])

        plt.plot(self.targ_mod.w, self.ref_mod.s - self.targ_mod.s,
                 '-', color='black')
        plots.annotate_spectrum(labels['residuals'], spec_offset=-1)

        ylim = plt.ylim(ymin=-0.5)
        minor_ticks = np.arange(ylim[0], ylim[1], 0.5)
        plt.yticks(minor_ticks)
        plt.grid(True, which='both')

    def to_hdf(self, outfile):
        """Saves the Match object to an hdf file.

        Args:
            outfile (str or h5 file): Output path or file handle.
        """
        # Allow either a string or h5 file object to be passed.
        is_path = False
        if isinstance(outfile, str):
            outfile = h5py.File(outfile, 'w')
            is_path = True

        # Specmatch-Emp version
        outfile['smemp-version'] = specmatchemp.SPECMATCH_VERSION.lstrip('v')

        # Save target
        outfile.create_group('target')
        self.target.to_hdf(outfile['target'])
        outfile.create_group('targ_mod')
        self.targ_mod.to_hdf(outfile['targ_mod'])

        # Save references
        outfile['num_refs'] = self.num_refs
        outfile.create_group('references')
        for i in range(self.num_refs):
            grp = outfile['references'].create_group('{0:d}'.format(i))
            self.refs[i].to_hdf(grp)

        outfile['ref_chisq'] = self.ref_chisq
        outfile['vsini'] = self.original_vsini

        # Save modified
        outfile.create_group('ref_mod')
        self.ref_mod.to_hdf(outfile['ref_mod'])

        # Save best-fit parameters
        outfile['best_params'] = self.best_params.dumps()
        outfile['spline'] = self.spl
        outfile['best_chisq'] = self.best_chisq

        if is_path:
            outfile.close()

    @classmethod
    def read_hdf(cls, infile):
        """Reads the Match object from an HDF file.

        Args:
            infile (str or h5 file): Input path or file handle.
        """
        # Allow either a string or h5 file object to be passed.
        is_path = False
        if isinstance(infile, str):
            infile = h5py.File(infile, 'r')
            is_path = True

        # Read target
        # Backward compatibility
        # if 'smemp-version' not in infile or infile['smemp-version'] < 0.4:
        #     # Read target
        #     target = spectrum.read_hdf(infile['target'])
        #     targ_mod = target.copy()
        #     ref_mod = spectrum.read_hdf(infile['modified'])
        # else:
        target = spectrum.read_hdf(infile['target'])
        targ_mod = spectrum.read_hdf(infile['targ_mod'])
        ref_mod = spectrum.read_hdf(infile['ref_mod'])

        # Read reference
        num_refs = infile['num_refs'].value
        ref_specs = []
        for i in range(num_refs):
            spec = spectrum.read_hdf(infile['references/{0:d}'.format(i)])
            ref_specs.append(spec)
        if 'ref_chisq' in infile:
            ref_chisq = infile['ref_chisq'][:]
        else:
            ref_chisq = None
        vsini = infile['vsini'][:]

        # Read best-fit parameters
        best_params = lmfit.Parameters()
        best_params.loads(infile['best_params'].value)
        spl = infile['spline'][:]
        best_chisq = infile['best_chisq'].value

        mt = cls(target, ref_specs, vsini)
        mt.load_params(best_params)
        mt.ref_chisq = ref_chisq

        if not np.allclose(ref_mod.s, mt.ref_mod.s) \
                or not np.allclose(targ_mod.s, mt.targ_mod.s) \
                or not np.allclose(spl, mt.spl) \
                or mt.best_chisq != best_chisq:
            print("Warning: Saved model and recreated model are not " +
                  "identical. Model may have been created by a different " +
                  "version of SpecMatch-Emp.")

        if is_path:
            infile.close()

        return mt


def add_spline_positions(params, knotx):
    """Adds spline positions to the parameter list.

    Args:
        params (lmfit.Parameters): parameters
        knotx (np.array): Array of knot positions
    Returns:
        params (lmfit.Parameters)
    """
    params.add('num_knots', value=len(knotx), vary=False)

    for i in range(len(knotx)):
        p = 'knotx_{0:d}'.format(i)
        params.add(p, value=knotx[i], vary=False)

    return params


def get_spline_positions(params):
    """Gets the spline positions from an lmfit parameters object.

    Args:
        params (lmfit.Parameters): parameters
    Returns:
        knotx (np.ndarray)
    """
    num_knots = params['num_knots'].value

    knotx = []
    for i in range(num_knots):
        p = 'knotx_{0:d}'.format(i)
        knotx.append(params[p].value)

    return np.array(knotx)


def add_vsini(params, vsini):
    """Adds vsini to an lmfit parameter list.

    Args:
        params (lmfit.Parameters): parameters
        vsini (np.ndarray): vsini values for each reference spectrum
    Returns:
        params (lmfit.Parameters)
    """
    if 'num_refs' not in params.valuesdict():
        params.add('num_refs', value=len(vsini), vary=False)

    for i in range(len(vsini)):
        p = 'vsini_{0:d}'.format(i)
        params.add(p, value=vsini[i], vary=False)

    return params


def get_vsini(params):
    """Gets vsini list from a parameters object.

    Args:
        params (lmfit.Parameters): parameters
    Returns:
        vsini (np.ndarray)
    """
    num_refs = params['num_refs'].value

    vsini = []
    for i in range(num_refs):
        p = 'vsini_{0:d}'.format(i)
        vsini.append(params[p].value)

    return np.array(vsini)


def add_lincomb_coeffs(params, num_refs):
    """Adds lincomb coefficients to an lmfit parameter list.

    Args:
        params (lmfit.Parameters): parameters
        num_refs (int): Number of reference spectra
    Returns:
        params (lmfit.Parameters)
    """
    if 'num_refs' not in params.valuesdict():
        params.add('num_refs', value=num_refs, vary=False)

    for i in range(num_refs):
        p = 'coeff_{0:d}'.format(i)
        params.add(p, value=1./num_refs, min=0.0, max=1.0)

    return params


def get_lincomb_coeffs(params):
    """Gets the lincomb coefficients form an lmfit parameter list.

    Args:
        params (lmfit.Paremters): parameters
    Returns:
        coeffs (np.ndarray)
    """
    num_refs = params['num_refs'].value

    coeffs = []
    for i in range(num_refs):
        p = 'coeff_{0:d}'.format(i)
        coeffs.append(params[p].value)

    return np.array(coeffs)


def get_physical_vsini(vsini):
    """Get the physical vsini

    Returns:
        physical vsini
    """
    if vsini >= 1.0:
        return vsini
    else:
        return -2.0 + vsini
