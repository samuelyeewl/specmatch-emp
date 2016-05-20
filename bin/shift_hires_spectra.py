#!/usr/bin/env python
"""
@filename shift_hires_spectra.py

Shift a target spectrum from HIRES onto a reference spectrum.
Refernce spectrum has to be in the standard SpecMatch format.
"""

from specmatchemp import shift_spectra
from specmatchemp import specmatch_io
import argparse
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Shift a target spectrum onto a reference spectrum')
    parser.add_argument('target_path', type=str)
    parser.add_argument('reference_path', type=str)
    parser.add_argument('output_path', type=str)
    args = parser.parse_args()

    s, w, serr = shift_spectra.main(args.target_path, 'hires', args.reference_path, args.output_path, True)
    
    nso_path = '/Users/samuel/Dropbox/SpecMatch-Emp/nso/nso_std.fits'
    s_nso, w_nso, serr_nso, h_nso = specmatch_io.read_standard_spectrum(nso_path)

    plt.clf()
    plt.plot(w_nso, s_nso)
    plt.plot(w, s)
    plt.plot(w, serr)
    plt.xlim(5180, 5200)
    plt.show()