#!/usr/bin/env python
"""
@filename plot_library.py

Plots global match plots
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

import re
import os
from argparse import ArgumentParser

from specmatchemp import library
from specmatchemp.plotting import plots

def main(libpath, outpath):
    lib = library.read_hdf(libpath, wavlim='none')

    with PdfPages(outpath) as pdf:
        fig = plt.figure(figsize=(12,8))
        plots.plot_library_params(lib, 'Teff', 'radius', grouped=True)
        plots.reverse_x()
        ax = plt.gca()
        ax.set_yscale('log')
        plt.ylim((0.1, 16))
        plt.legend(loc='upper left')
        pdf.savefig()
        plt.close()

        fig = plt.figure(figsize=(12,8))
        plots.plot_library_params(lib, 'feh', 'radius', grouped=True)
        ax = plt.gca()
        ax.set_yscale('log')
        plt.ylim((0.1, 16))
        plt.legend(loc='upper left')
        pdf.savefig()
        plt.close()


if __name__ == '__main__':
    psr = ArgumentParser(description="Produce plots for analysis of match results")
    psr.add_argument('library', type=str, help="Path to library h5 file")
    psr.add_argument('outpath', type=str, help="Path to output file")
    args = psr.parse_args()

    main(args.library, args.outpath)