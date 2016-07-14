#!/usr/bin/env python

import pandas as pd
import glob
import os
import sys
from argparse import ArgumentParser

if __name__ == '__main__':
    # Argument parser
    psr = ArgumentParser(description="Combine match results")
    psr.add_argument('dir', type=str, help="Path to directory")
    psr.add_argument('suffix', type=str, help="File suffix to search for")
    psr.add_argument('outfile', type=str, help="Output filename")
    args = psr.parse_args()

    if not os.path.isdir(args.dir):
        print("Could not find folder at {0}".format(args.dir))
        sys.exit()

    cols = ['targ_idx', 'ref_idx', 'chi_squared', 'fit_params']
    df = pd.DataFrame(columns=cols)

    for file in glob.glob(os.path.join(args.dir, "*"+ args.suffix)):
        res = pd.DataFrame.from_csv(file)
        df = pd.concat((df, res))

    df.reset_index(inplace=True)
    df.to_csv(args.outfile)
