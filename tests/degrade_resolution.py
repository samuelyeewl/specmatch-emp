#!/usr/bin/env python
"""
@filename degrade_resolution.py

Degrade each order of a HIRES spectrum to the given resolution.
"""
import numpy as np
import os
from argparse import ArgumentParser

from scipy.signal import gaussian
from scipy.ndimage.filters import convolve1d

from specmatchemp import SPECMATCHDIR
from specmatchemp import spectrum


def degrade_spec(spec, r):
    """Degrade the spectrum to the given spectral resolution by convolving with
    a Gaussian kernel.

    Args:
        spec (spectrum.HiresSpectrum): Target spectrum
        r (int): Desired spectral resolution
    """
    spec = spec.copy()

    num_orders = spec.s.shape[0]
    for i in range(num_orders):
        # Perform convolution on every order
        mean_w = np.mean(spec.w[i])
        fwhm = mean_w / r
        # Std dev of Gaussian = FWHM / 2.355
        sigma_w = fwhm / 2.355
        # Width of Gaussian in pixels
        dw = np.mean(spec.w[i, 1:] - spec.w[i, :-1])
        sigma_pix = sigma_w / dw

        # number of pixels
        n = 201
        # Create kernel
        kernel = gaussian(n, sigma_pix)
        kernel /= np.sum(kernel)

        # perform convolution
        spec.s[i] = convolve1d(spec.s[i], kernel)
        spec.serr[i] = convolve1d(spec.serr[i], kernel)

    return spec


def main(args):
    filename = os.path.splitext(args.input_file)[0]
    infile = os.path.join(args.dir, filename + '.fits')

    spec = spectrum.read_hires_fits(infile)

    d_spec = degrade_spec(spec, args.R)

    suffix = "_R={0:d}".format(args.R)
    outfile = os.path.join(args.dir, filename + suffix + '.fits')

    # Save degraded spectrum
    d_spec.to_hires_fits(outfile, clobber=True)


if __name__ == '__main__':
    # Argument parser
    psr = ArgumentParser(description="Degrade the resolution of a given spectrum")
    psr.add_argument('input_file', type=str, help="Filename of input spectrum")
    psr.add_argument('R', type=int, help="Target spectral resolution")
    psr.add_argument('-d', '--dir', type=str, default=SPECMATCHDIR+'spectra/',
                     help="Directory to look for spectra")
    args = psr.parse_args()

    main(args)
