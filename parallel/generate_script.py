#!/usr/bin/env python
"""
@filename generate_plot_script.py

Generate script lines for shifting spectra
"""

import os.path, sys
import pandas as pd
from argparse import ArgumentParser
from specmatchemp import library

if __name__ == '__main__':
    psr = ArgumentParser(description="Generic script generation")
    psr.add_argument('outpath', type=str, help="Path to save script")
    psr.add_argument('libpath', type=str, help="Path to library file")
    psr.add_argument('formatstring', type=str, help="Format string")
    psr.add_argument('cols', nargs='*', type=str, help="Columns to use")
    psr.add_argument('-a', '--append', action='store_const', const='a+', default="w+", help="Append to existing script")
    args = psr.parse_args()

    if not os.path.isfile(args.libpath):
        print("Could not find {0}".format(args.libpath))
        sys.exit(1)

    if os.path.splitext(args.libpath)[1] == ".h5":
        lib = library.read_hdf(args.libpath, wavlim='none')
        df = lib.library_params
    elif os.path.splitext(args.libpath)[1] == ".csv":
        df = pd.read_csv(args.libpath)
    else:
        print("{0} was not a valid library format".format(args.libpath))

    scriptfile = open(args.outpath, args.append)

    for idx, row in df.iterrows():
        s = "source ~/.bash_profile; "
        s+= args.formatstring.format(*(row[args.cols].values))
        s+= "\n"
        scriptfile.write(s)

    scriptfile.close()






