#!/usr/bin/env python
"""
@filename get_isochrones.py

Builds the parameter table of the SpecMatch-Emp library by reading in the various catalogs, checking
for spectra in the CPS database.
Uses the isochrones package to obtain a full set of stellar parameters from those which are known.
Saves the stars and parameters as a Pandas Dataframe.
"""

from __future__ import print_function

import sys, os
from argparse import ArgumentParser

import numpy as np
import pandas as pd

from isochrones.dartmouth import Dartmouth_Isochrone
from isochrones import StarModel

from specmatchemp import library


def get_isochrone_params(stars, modeldir, overwrite=False):
    """Fill out parameter table with values obtained from isochrone package

    Args:
        stars (pd.DataFrame): parameter table
        modeldir (str)   : directory to save fitted models
        overwrite (bool) : (optional) whether to use existing models or overwrite
    Returns:
        stars (pd.DataFrame): updated parameter table
        bad_idxs (array)    : indices of stars with multiple parameters inconsistent with model
    """
    dar = Dartmouth_Isochrone()

    num_stars = len(stars)
    current_star = 1
    bad_idxs = []
    for i, row in stars.iterrows():
        print("Getting isochrone parameters for star {0} of {1}".format(current_star, num_stars))
        current_star += 1

        # get known stellar properties
        lib_props = {}
        for p in library.STAR_PROPS:
            if not np.isnan(row[p]):
                # (value, uncertainty)
                lib_props[p] = (row[p], row['u_'+p])
        # if all properties are known, we don't need to use the model
        if len(lib_props) == 6:
            continue

        # isochrones requires fitting a distance to at least one magnitude - use arbitrary value
        lib_props['V'] = 1.

        if modeldir is None:
            model = StarModel(dar, **lib_props)
            model.fit(overwrite=True, verbose=False)
        else:
            # check if fitting has already been done
            modelfile = os.path.join(modeldir, "{0}_model.h5".format(row['cps_name']))
            if os.path.exists(modelfile) and not overwrite:
                model = StarModel.load_hdf(modelfile)
            # otherwise perform the fit
            else:
                print(len(lib_props))
                print(lib_props)
                print(row)
                raise Exception
                model = StarModel(dar, **lib_props)
                model.fit(overwrite=True, verbose=False)
                model.save_hdf(modelfile)

        N_SIGMA = 2
        MAX_PERCENTILE = 0.95
        MIN_PERCENTILE = 0.05
        incons_props = 0
        # fill out unknown parameters
        for p in library.STAR_PROPS:
            value = model.samples[p].quantile(0.5)
            upper_bound = model.samples[p].quantile(MAX_PERCENTILE)
            lower_bound = model.samples[p].quantile(MIN_PERCENTILE)
            if p in lib_props:
                # check for model consistency with known library values
                # check if 2-sigma bounds overlap
                if (row[p]+N_SIGMA*row['u_'+p]) < lower_bound or (row[p]-N_SIGMA*row['u_'+p]) > upper_bound:
                    print("Warning: Model for star {0} had inconsistent values in {1}.".format(
                        row['cps_name'], p))
                    # Save error messages to file
                    outpath_err = os.path.join(modeldir, 'model_inconsistencies.txt')
                    with open(outpath_err, 'a') as f:
                        f.write("Inconsistent {0} for star {1}\n".format(p, row['cps_name']))
                        f.write("\tLibrary values: {0:.2f} +/- {1:.2f}\n".format(row[p], row['u_'+p]))
                        f.write("\tModel values: {0:.2f}, {1:d}-sigma = ({2:.2f}, {3:.2f})\n".format(
                        value, N_SIGMA, lower_bound, upper_bound))
                    incons_props += 1
            else:
                # insert unknown values
                stars.loc[i, p] = np.around(value, 2)
                stars.loc[i, 'u_'+p] = np.around((upper_bound-lower_bound)/2, 2)

        # record index if more than 2 properties were inconsistent with library values
        if incons_props > 1:
            bad_idxs.append(i)

    return stars, bad_idxs


def main(libpath, modeldir, overwrite, maskpath=''):
    ext = os.path.splitext(libpath)[-1].lower()
    if ext == '.csv':
        stars = pd.read_csv(libpath, index_col=0)
        mask = pd.read_csv(maskpath, index_col=0)
        stars, bad_idxs = get_isochrone_params(stars, modeldir, overwrite)
        # drop bad indexes
        stars.drop(bad_idxs, inplace=True)
        mask.drop(bad_idxs, inplace=True)
        stars.to_csv(libpath)
        mask.to_csv(maskpath)
    elif ext == '.h5':
        lib = library.read_hdf(libpath)
        lib.library_params, bad_idxs = get_isochrone_params(lib.library_params, modeldir, overwrite)
        # drop bad indexes
        bad_idxs.sort(reverse=True)
        for idx in bad_idxs:
            lib.remove(idx)
        lib.to_hdf(libpath)
    else:
        print("Library file was not in a valid format.")
        sys.exit(1)


if __name__ == '__main__':
    psr = ArgumentParser(description="Fill in unknown parameters by using isochrone models")
    psr.add_argument('libpath', type=str, help="Path to dataframe (.csv) or specmatchemp.library file (.h5)")
    psr.add_argument('modeldir', type=str, nargs='?', default=None, help="Directory to store models in")
    psr.add_argument('-o', '--overwrite', action='store_true', help='Generate models even if they exist')
    psr.add_argument('-m', '--maskpath', type=str, help='Path to mask dataframe (.csv)')
    args = psr.parse_args()

    main(args.libpath, args.modeldir, args.overwrite, args.maskpath)
