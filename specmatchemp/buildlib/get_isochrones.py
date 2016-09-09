#!/usr/bin/env python
"""
@filename get_isochrones.py

Uses the isochrones package to obtain a full set of stellar parameters from
those which are known.
Also creates a mask indicating which parameters are from the literature (True)
and which are obtained using isochrones.
"""

from __future__ import print_function

import os
import sys
from argparse import ArgumentParser

import numpy as np
import pandas as pd

from isochrones.dartmouth import Dartmouth_Isochrone
from isochrones import StarModel

from specmatchemp import SPECMATCHDIR
from specmatchemp import library


def get_isochrone_params(stars, modeldir, overwrite=False):
    """Fill out parameter table with values obtained from isochrone package

    Args:
        stars (pd.DataFrame): parameter table
        modeldir (str): directory to save fitted models
        overwrite (bool, optional): whether to use existing models or overwrite
    Returns:
        stars (pd.DataFrame): updated parameter table
    """
    dar = Dartmouth_Isochrone()

    num_stars = len(stars)
    current_star = 1

    for idx, row in stars.iterrows():
        print("Getting isochrone parameters for star {0} of {1}"
              .format(current_star, num_stars))
        current_star += 1

        # get known stellar properties
        lib_props = {}
        for p in library.Library.STAR_PROPS:
            if not np.isnan(row[p]):
                # (value, uncertainty)
                lib_props[p] = (row[p], row['u_'+p])
        # if all properties are known, we don't need to use the model
        if len(lib_props) == 6:
            continue

        if modeldir is None:
            model = StarModel(dar, **lib_props)
            model.fit(overwrite=True, verbose=False)
        else:
            # check if fitting has already been done
            modelfile = os.path.join(modeldir,
                                     "{0}_model.h5".format(row['cps_name']))
            if os.path.exists(modelfile) and not overwrite:
                model = StarModel.load_hdf(modelfile)
            # otherwise perform the fit
            else:
                model = StarModel(dar, **lib_props)
                model.fit(overwrite=True, verbose=False)
                model.save_hdf(modelfile)

        N_SIGMA = 2
        MAX_PERCENTILE = 0.95
        MIN_PERCENTILE = 0.05
        # fill out unknown parameters
        for p in library.Library.STAR_PROPS:
            value = model.samples[p].quantile(0.5)
            upper_bound = model.samples[p].quantile(MAX_PERCENTILE)
            lower_bound = model.samples[p].quantile(MIN_PERCENTILE)

            # If a property is already known, check for model consistency
            if p in lib_props:
                # check if 2-sigma bounds fail to overlap
                if (row[p] + N_SIGMA * row['u_'+p]) < lower_bound \
                        or (row[p] - N_SIGMA * row['u_'+p]) > upper_bound:
                    warningstr = "Warning: Inconisistent {0} for star {1}\n"\
                        .format(p, row['cps_name'])
                    warningstr += "\tLibrary values: {0:.2f} +/- {1:.2f}\n"\
                        .format(row[p], row['u_'+p])
                    warningstr += "\tModel values: {0:.2f}, ".format(value)
                    warningstr += "{0:d}-sigma = ({1:.2f}, {2:.2f})\n".format(
                        N_SIGMA, lower_bound, upper_bound)

                    print(warningstr)
                    # Save error messages to file
                    errpath = os.path.join(modeldir, 'errors.txt')
                    with open(errpath, 'a') as f:
                        f.write(warningstr)
            # Insert the unknown values if we don't already know them
            else:
                stars.loc[idx, p] = np.around(value, 2)
                stars.loc[idx, 'u_'+p] = \
                    np.around((upper_bound - lower_bound) / 2, 2)
    return stars


def main(libpath, modeldir, overwrite):
    ext = os.path.splitext(libpath)[-1].lower()
    if ext == '.csv':
        params = pd.read_csv(libpath)
        maskpath = libpath[:-4] + '_mask.csv'
        # Create mask dataframe
        param_mask = pd.DataFrame()
        for p in library.Library.STAR_PROPS:
            param_mask[p] = np.isfinite(params[p])
            param_mask['u_'+p] = np.isfinite(params[p])

        params = get_isochrone_params(params, modeldir, overwrite)

        params.to_csv(libpath)
        param_mask.to_csv(maskpath)

    elif ext == '.h5':
        lib = library.read_hdf(libpath)
        param_mask = pd.DataFrame()
        for p in library.Library.STAR_PROPS:
            param_mask[p] = np.isfinite(lib.library_params[p])
            param_mask['u_'+p] = np.isfinite(lib.library_params[p])

        lib.library_params = get_isochrone_params(lib.library_params,
                                modeldir, overwrite)
        lib.param_mask = param_mask

        lib.to_hdf(libpath)

    else:
        print("Library file was not in a valid format.")
        sys.exit(1)


if __name__ == '__main__':
    psr = ArgumentParser(description="Fill in unknown parameters by using " +
                                     "isochrone models")
    psr.add_argument('-l', '--libpath', type=str,
                     default=os.path.join(SPECMATCHDIR, 'libstars.csv'),
                     help="Path to dataframe (.csv) or specmatchemp.library " +
                     "file (.h5)")
    psr.add_argument('-m', '--modeldir', type=str,
                     default=os.path.join(SPECMATCHDIR, 'isochrone_models/'),
                     help="Directory to store isochrone models")
    psr.add_argument('-o', '--overwrite', action='store_true',
                     help='Generate models even if they exist')
    args = psr.parse_args()

    main(args.libpath, args.modeldir, args.overwrite)
