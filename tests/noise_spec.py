#!/usr/bin/env python
"""
@filename noise_spec.py

Adds noise to a given spectrum
"""
import numpy as np
import os
from argparse import ArgumentParser

from specmatchemp import SPECMATCHDIR
from specmatchemp import spectrum


def add_noise(spec, snr):
    """Adds Poisson noise to the given spectrum.

    The amount of noise added is such that when added in quadrature with the
    current noise level, the target signal-to-noise ratio is achieved.

    Args:
        spec (spectrum.Spectrum): Target spectrum
        snr (float): Target SNR to achieve
    """
    spec = spec.copy()
    current_snr = spec.snr()
    extra_error = np.sqrt((1/snr)**2 - (1/current_snr)**2)

    # Use Gaussian approximation to generate noise
    noise = np.random.normal(scale=extra_error, size=spec.s.shape)

    # Add noise to spectrum, scalled accordingly
    spec.s += spec.s * noise

    return spec


def main(args):
    filename = os.path.splitext(args.input_file)[0]
    infile = os.path.join(args.dir, filename + '.fits')

    spec = spectrum.read_hires_fits(infile)

    for i in range(args.num):
        noised_spec = add_noise(spec, args.snr)

        if args.num == 1:
            suffix = "_snr={0:.0f}".format(args.snr)
        else:
            suffix = "_snr={0:.0f}_i={1:d}".format(args.snr, i + 1)
        outfile = os.path.join(args.dir, filename + suffix + '.fits')

        noised_spec.to_hires_fits(outfile)


if __name__ == '__main__':
    # Argument parser
    psr = ArgumentParser(description="Add Poisson noise to a spectrum")
    psr.add_argument('input_file', type=str, help="Filename of input spectrum")
    psr.add_argument('snr', type=float, help="Target signal-to-noise ratio")
    psr.add_argument('-n', '--num', type=int, default=1,
                     help="Number of noisy spectra to generate")
    psr.add_argument('-d', '--dir', type=str, default=SPECMATCHDIR+'spectra/',
                     help="Directory to look for spectra")
    args = psr.parse_args()

    main(args)
