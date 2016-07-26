#!/usr/bin/env python
"""
@filename read_catalogs.py

Builds the parameter table of the SpecMatch-Emp library by reading in the various catalogs, checking
for spectra in the CPS database.
Uses the isochrones package to obtain a full set of stellar parameters from those which are known.
Saves the stars and parameters as a Pandas Dataframe.
"""

import os
from argparse import ArgumentParser

import numpy as np
import pandas as pd

from astropy.io import ascii
from isochrones.dartmouth import Dartmouth_Isochrone
from isochrones import StarModel

from specmatchemp.library import LIB_COLS
from specmatchemp.library import STAR_PROPS
from specmatchemp.utils import cpsutils
from specmatchemp.utils import utils

# relative catalog locations
MANN_FILENAME = "Mann2015/stars.dat"
MANN_README = "Mann2015/ReadMe"
HUBER_FILENAME = "Huber2013/table2.dat"
HUBER_README = "Huber2013/ReadMe"
VONBRAUN_FILENAME = "VonBraun-Figure6.txt"
BREWER_FILENAME = "brewer_cut.csv"
RAMIREZ_FILENAME = "ramirez_2005_cut.csv"
CASAGRANDE_FILENAME = "Casagrande2006/table1_cut.csv"
BRUNTT_FILENAME = "Bruntt2012/table3.dat"
BRUNTT_README = "Bruntt2012/ReadMe.txt"
CPS_INDEX = "cps_templates.csv"
CPS_SPECTRA_DIR = "iodfitsdb/"

NOSPECTRA_COLS = ['name', 'source']

MIN_PERCENTILE = 0.05
MAX_PERCENTILE = 0.95

def read_brewer(catalogdir, cps_list):
    """Read in Brewer (2016) catalog

    Args:
        catalogdir (str): Path of catalog directory.
        cps_list (pd.DataFrame): Dataframe containing list of CPS spectra.

    Returns:
        stars (pd.DataFrame): Stars in source which have CPS spectra
        nospec (pd.DataFrame): Stars in source which don't have CPS spectra
    """
    brewer_data = ascii.read(os.path.join(catalogdir, BREWER_FILENAME))

    stars = pd.DataFrame(columns=LIB_COLS)
    nospectra = pd.DataFrame(columns=NOSPECTRA_COLS)

    for row in brewer_data:
        try:
            query_result = cpsutils.find_spectra(row['NAME'], cps_list)
            if not query_result.empty:
                new_row = {}
                new_row['cps_name'] = str(query_result.iloc[0]['name'])
                new_row['obs'] = query_result.obs.values
                new_row['Teff'] = row['TEFF']
                new_row['u_Teff'] = 25
                new_row['logg'] = row['LOGG']
                new_row['u_logg'] = 0.028
                new_row['feh'] = row['FEH']
                new_row['u_feh'] = 0.010
                new_row['vsini'] = row['VSINI']
                new_row['source'] = 'Brewer'
                new_row['source_name'] = row['NAME']

                stars = stars.append(pd.Series(new_row), ignore_index=True)
            else:
                new_row = {}
                new_row['name'] = row['NAME']
                new_row['source'] = 'Brewer'
                nospectra = nospectra.append(pd.Series(new_row), ignore_index=True)
        except:
            new_row = {}
            new_row['name'] = row['NAME']
            new_row['source'] = 'Brewer'
            nospectra = nospectra.append(pd.Series(new_row), ignore_index=True)

    return stars, nospectra

def read_mann(catalogdir, cps_list):
    """Read in Mann (2015) catalog

    Args:
        catalogdir (str): Path of catalog directory.
        cps_list (pd.DataFrame): Dataframe containing list of CPS spectra.

    Returns:
        stars (pd.DataFrame): Stars in source which have CPS spectra
        nospec (pd.DataFrame): Stars in source which don't have CPS spectra
    """
    mann_data = ascii.read(os.path.join(catalogdir, MANN_FILENAME), readme=os.path.join(catalogdir, MANN_README))

    stars = pd.DataFrame(columns=LIB_COLS)
    nospectra = pd.DataFrame(columns=NOSPECTRA_COLS)

    for row in mann_data:
        try:
            query_result = cpsutils.find_spectra(row['CNS3'], cps_list)
            if not query_result.empty:
                new_row = {}
                new_row['cps_name'] = str(query_result.iloc[0]['name'])
                new_row['obs'] = query_result.obs.values
                new_row['Teff'] = row['Teff']
                new_row['u_Teff'] = row['e_Teff']
                new_row['radius'] = row['R']
                new_row['u_radius'] = row['e_R']
                new_row['feh'] = row['[Fe/H]']
                new_row['u_feh'] = row['e_[Fe/H]']
                new_row['mass'] = row['M']
                new_row['u_mass'] = row['e_M']
                new_row['source'] = 'Mann'
                new_row['source_name'] = row['CNS3']
                # Calculate logg from radius and mass
                logg, u_logg = utils.calc_logg(row['R'], row['e_R'], row['M'], row['e_M'])
                new_row['logg'] = logg
                new_row['u_logg'] = u_logg

                stars = stars.append(pd.Series(new_row), ignore_index=True)
            else:
                new_row = {}
                new_row['name'] = row['CNS3']
                new_row['source'] = 'Mann'
                nospectra = nospectra.append(pd.Series(new_row), ignore_index=True)
        except:
            new_row = {}
            new_row['name'] = row['CNS3']
            new_row['source'] = 'Mann'
            nospectra = nospectra.append(pd.Series(new_row), ignore_index=True)

    return stars, nospectra

def read_vonbraun(catalogdir, cps_list):
    """Read in Von Braun (2013) catalog

    Args:
        catalogdir (str): Path of catalog directory.
        cps_list (pd.DataFrame): Dataframe containing list of CPS spectra.

    Returns:
        stars (pd.DataFrame): Stars in source which have CPS spectra
        nospec (pd.DataFrame): Stars in source which don't have CPS spectra
    """
    vb_data = ascii.read(os.path.join(catalogdir, VONBRAUN_FILENAME))

    stars = pd.DataFrame(columns=LIB_COLS)
    nospectra = pd.DataFrame(columns=NOSPECTRA_COLS)

    for row in vb_data:
        try:
            query_result = cpsutils.find_spectra(row['Star'], cps_list)
            if not query_result.empty:
                new_row = {}
                new_row['cps_name'] = str(query_result.iloc[0]['name'])
                new_row['obs'] = query_result.obs.values
                new_row['Teff'] = row['Teff']
                new_row['u_Teff'] = row['eTeff']
                new_row['radius'] = row['Radius']
                new_row['u_radius'] = row['eRadius']
                new_row['feh'] = row['FeH']
                new_row['u_feh'] = 0.10
                new_row['source'] = 'Von Braun'
                new_row['source_name'] = row['Star']

                stars = stars.append(pd.Series(new_row), ignore_index=True)
            else:
                new_row = {}
                new_row['name'] = row['Star']
                new_row['source'] = 'Von Braun'
                nospectra = nospectra.append(pd.Series(new_row), ignore_index=True)
        except:
            new_row = {}
            new_row['name'] = row['Star']
            new_row['source'] = 'Von Braun'
            nospectra = nospectra.append(pd.Series(new_row), ignore_index=True)

    return stars, nospectra

def read_huber(catalogdir, cps_list):
    """Read in Huber (2013) catalog

    Args:
        catalogdir (str): Path of catalog directory.
        cps_list (pd.DataFrame): Dataframe containing list of CPS spectra.

    Returns:
        stars (pd.DataFrame): Stars in source which have CPS spectra
        nospec (pd.DataFrame): Stars in source which don't have CPS spectra
    """
    huber_data = ascii.read(os.path.join(catalogdir, HUBER_FILENAME), readme=os.path.join(catalogdir, HUBER_README))
    huber_data = huber_data[huber_data['f_KOI'] != '*']

    stars = pd.DataFrame(columns=LIB_COLS)
    nospectra = pd.DataFrame(columns=NOSPECTRA_COLS)

    for row in huber_data:
        try:
            query_result = cpsutils.find_spectra('KOI'+str(row['KOI']), cps_list)
            if not query_result.empty:
                new_row = {}
                new_row['cps_name'] = str(query_result.iloc[0]['name'])
                new_row['obs'] = query_result.obs.values
                new_row['Teff'] = row['Teff']
                new_row['u_Teff'] = row['e_Teff']
                new_row['radius'] = float(row['Rad'])
                new_row['u_radius'] = row['e_Rad']
                new_row['feh'] = row['[Fe/H]']
                new_row['u_feh'] = row['e_[Fe/H]']
                new_row['mass'] = row['Mass']
                new_row['u_mass'] = row['e_Mass']
                new_row['source'] = 'Huber'
                new_row['source_name'] = 'KOI'+str(row['KOI'])
                # Calculate logg from radius and mass
                logg, u_logg = utils.calc_logg(row['Rad'], row['e_Rad'], row['Mass'], row['e_Mass'])
                new_row['logg'] = logg
                new_row['u_logg'] = u_logg

                stars = stars.append(pd.Series(new_row), ignore_index=True)
            else:
                new_row = {}
                new_row['name'] = 'KOI'+str(row['KOI'])
                new_row['source'] = 'Huber'
                nospectra = nospectra.append(pd.Series(new_row), ignore_index=True)
        except:
            new_row = {}
            new_row['name'] = 'KOI'+str(row['KOI'])
            new_row['source'] = 'Huber'
            nospectra = nospectra.append(pd.Series(new_row), ignore_index=True)

    return stars, nospectra

def read_ramirez(catalogdir, cps_list):
    """Read in Ramirez (2005) catalog

    Args:
        catalogdir (str): Path of catalog directory.
        cps_list (pd.DataFrame): Dataframe containing list of CPS spectra.

    Returns:
        stars (pd.DataFrame): Stars in source which have CPS spectra
        nospec (pd.DataFrame): Stars in source which don't have CPS spectra
    """
    ramirez_data = pd.read_csv(os.path.join(catalogdir, RAMIREZ_FILENAME))

    stars = pd.DataFrame(columns=LIB_COLS)
    nospectra = pd.DataFrame(columns=NOSPECTRA_COLS)

    for idx, row in ramirez_data.iterrows():
        try:
            query_result = cpsutils.find_spectra(row['name'], cps_list)
            if not query_result.empty:
                new_row = {}
                new_row['cps_name'] = str(query_result.iloc[0]['name'])
                new_row['obs'] = query_result.obs.values
                new_row['Teff'] = row['Teff']
                new_row['u_Teff'] = row['u_Teff']
                new_row['logg'] = row['logg']
                new_row['u_logg'] = row['u_logg']
                new_row['feh'] = row['feh']
                new_row['u_feh'] = row['u_feh']
                new_row['source'] = 'Ramirez'
                new_row['source_name'] = row['name']

                stars = stars.append(pd.Series(new_row), ignore_index=True)
            else:
                new_row = {}
                new_row['name'] = row['name']
                new_row['source'] = 'Ramirez'
                nospectra = nospectra.append(pd.Series(new_row), ignore_index=True)
        except:
            new_row = {}
            new_row['name'] = row['name']
            new_row['source'] = 'Ramirez'
            nospectra = nospectra.append(pd.Series(new_row), ignore_index=True)

    return stars, nospectra

def read_casagrande(catalogdir, cps_list):
    """Read in Casagrande (2006) catalog

    Args:
        catalogdir (str): Path of catalog directory.
        cps_list (pd.DataFrame): Dataframe containing list of CPS spectra.

    Returns:
        stars (pd.DataFrame): Stars in source which have CPS spectra
        nospec (pd.DataFrame): Stars in source which don't have CPS spectra
    """
    c_data = pd.read_csv(os.path.join(catalogdir, CASAGRANDE_FILENAME))

    stars = pd.DataFrame(columns=LIB_COLS)
    nospectra = pd.DataFrame(columns=NOSPECTRA_COLS)

    for idx, row in c_data.iterrows():
        try:
            query_result = cpsutils.find_spectra(row['Name'], cps_list)
            if not query_result.empty:
                # Calculate stellar radius from angular diameter and parallax
                radius, u_radius = utils.calc_radius(row['Plx'], row['e_Plx'], row['Diam'], row['e_Diam'])
                # Exclude stars with u_Teff > 100 or u_R/R > 0.07
                if row['e_Teff'] > 100 or u_radius/radius > 0.07:
                    continue
                new_row = {}
                new_row['cps_name'] = str(query_result.iloc[0]['name'])
                new_row['obs'] = query_result.obs.values
                new_row['Teff'] = row['Teff']
                new_row['u_Teff'] = row['e_Teff']
                new_row['feh'] = row['[Fe/H]']
                new_row['u_feh'] = 0.15
                new_row['radius'] = radius
                new_row['u_radius'] = u_radius
                new_row['source'] = 'Casagrande'
                new_row['source_name'] = row['Name']
                stars = stars.append(pd.Series(new_row), ignore_index=True)
            else:
                new_row = {}
                new_row['name'] = row['Name']
                new_row['source'] = 'Casagrande'
                nospectra = nospectra.append(pd.Series(new_row), ignore_index=True)
        except:
            new_row = {}
            new_row['name'] = row['name']
            new_row['source'] = 'Casagrande'
            nospectra = nospectra.append(pd.Series(new_row), ignore_index=True)

    return stars, nospectra

def read_bruntt(catalogdir, cps_list):
    """Read in Bruntt (2012) catalog

    Args:
        catalogdir (str): Path of catalog directory.
        cps_list (pd.DataFrame): Dataframe containing list of CPS spectra.

    Returns:
        stars (pd.DataFrame): Stars in source which have CPS spectra
        nospec (pd.DataFrame): Stars in source which don't have CPS spectra
    """
    b_data = ascii.read(os.path.join(catalogdir, BRUNTT_FILENAME), readme=os.path.join(catalogdir, BRUNTT_README))

    stars = pd.DataFrame(columns=LIB_COLS)
    nospectra = pd.DataFrame(columns=NOSPECTRA_COLS)

    for row in b_data:
        try:
            query_result = cpsutils.find_spectra('KIC'+str(row['KIC']), cps_list)
            if not query_result.empty:
                new_row = {}
                new_row['cps_name'] = str(query_result.iloc[0]['name'])
                new_row['obs'] = query_result.obs.values
                new_row['Teff'] = row['Teff']
                new_row['u_Teff'] = 60
                new_row['feh'] = row['[Fe/H]']
                new_row['u_feh'] = 0.06
                new_row['logg'] = row['logg']
                new_row['u_logg'] = 0.03
                new_row['vsini'] = row['vsini']
                new_row['source'] = 'Bruntt'
                new_row['source_name'] = 'KIC'+str(row['KIC'])

                stars = stars.append(pd.Series(new_row), ignore_index=True)
            else:
                new_row = {}
                new_row['name'] = 'KIC'+str(row['KIC'])
                new_row['source'] = 'Bruntt'
                nospectra = nospectra.append(pd.Series(new_row), ignore_index=True)
        except:
            new_row = {}
            new_row['name'] = 'KIC'+str(row['KIC'])
            new_row['source'] = 'Bruntt'
            nospectra = nospectra.append(pd.Series(new_row), ignore_index=True)


    return stars, nospectra

def read_catalogs(catalogdir, cpsdir):
    """Reads in the catalogs

    Args:
        catalogdir (str): Path to catalog directory
        cpsdir (str): Path to CPS list

    Returns:
        stars (pd.DataFrame): Stars in source which have CPS spectra
        nospec (pd.DataFrame): Stars in source which don't have CPS spectra
    """

    # Read in list of CPS spectra
    cps_index_path = os.path.join(cpsdir, CPS_INDEX)
    cps_list = pd.read_csv(cps_index_path)

    # Create dataframe to store found stars
    stars = pd.DataFrame(columns=LIB_COLS)

    # Create dataframe to store stars which have no spectra
    stars_nospectra = pd.DataFrame(columns=NOSPECTRA_COLS)

    # Read catalogs
    print("Reading Brewer catalog")
    brewer_stars, brewer_nospec = read_brewer(catalogdir, cps_list)
    stars = pd.concat((stars, brewer_stars), ignore_index=True)
    stars_nospectra = stars_nospectra.append(brewer_nospec)
    print("\t{0:d} stars with spectra from Brewer catalog".format(len(brewer_stars)))

    print("Reading Mann catalog")
    mann_stars, mann_nospec = read_mann(catalogdir, cps_list)
    stars = pd.concat((stars, mann_stars), ignore_index=True)
    stars_nospectra = stars_nospectra.append(mann_nospec)
    print("\t{0:d} stars with spectra from Mann catalog".format(len(mann_stars)))

    print("Reading von Braun catalog")
    vonbraun_stars, vonbraun_nospec = read_vonbraun(catalogdir, cps_list)
    stars = pd.concat((stars, vonbraun_stars), ignore_index=True)
    stars_nospectra = stars_nospectra.append(vonbraun_nospec)
    print("\t{0:d} stars with spectra from von Braun catalog".format(len(vonbraun_stars)))

    # print("Reading Huber catalog")
    # huber_stars, huber_nospec = read_huber(catalogdir, cps_list)
    # stars = pd.concat((stars, huber_stars), ignore_index=True)
    # stars_nospectra = stars_nospectra.append(huber_nospec)
    # print("\t{0:d} stars with spectra from Huber catalog".format(len(huber_stars)))

    print("Reading Ramirez catalog")
    ramirez_stars, ramirez_nospec = read_ramirez(catalogdir, cps_list)
    stars = pd.concat((stars, ramirez_stars), ignore_index=True)
    stars_nospectra = stars_nospectra.append(ramirez_nospec)
    print("\t{0:d} stars with spectra from Ramirez catalog".format(len(ramirez_stars)))

    print("Reading Casagrande catalog")
    c_stars, c_nospec = read_casagrande(catalogdir, cps_list)
    stars = pd.concat((stars, c_stars), ignore_index=True)
    stars_nospectra = stars_nospectra.append(c_nospec)
    print("\t{0:d} stars with spectra from Casagrande catalog".format(len(c_stars)))

    print("Reading Bruntt catalog")
    b_stars, b_nospec = read_bruntt(catalogdir, cps_list)
    stars = pd.concat((stars, b_stars), ignore_index=True)
    stars_nospectra = stars_nospectra.append(b_nospec)
    print("\t{0:d} stars with spectra from Bruntt catalog".format(len(b_stars)))

    dups = stars[stars.duplicated(subset='cps_name', keep=False)].sort_values(by='cps_name')
    dups_vb = dups[~dups.source.str.contains('Von Braun')]
    idxs = dups_vb.index
    print("Removing {0:d} duplicates, favoring von Braun data".format(len(idxs)))
    stars.drop(idxs, inplace=True)

    print("Total of {0:d} stars read".format(len(stars)))

    return stars, stars_nospectra


def get_isochrone_params(stars, diagnostic=False, outdir='./'):
    """Fill out parameter table with values obtained from isochrone package

    Args:
        stars (pd.DataFrame): star library
        diagnostic (bool)   : whether to save the fitted model as a file
        outdir (str)        : directory to save fitted models
    Returns:
        stars (pd.DataFrame): star library with updated parameters
    """
    dar = Dartmouth_Isochrone()

    # create output subfolder
    if diagnostic:
        outdir = os.path.join(outdir, 'isochrone_models/')
        if not os.path.exists(outdir):
            os.makedirs(outdir)

    num_stars = len(stars)
    current_star = 1

    for i, row in stars.iterrows():
        print("Getting isochrone parameters for star {0} of {1}".format(current_star, num_stars))
        current_star += 1
        # get known stellar properties
        lib_props = {}
        for p in STAR_PROPS:
            if not np.isnan(row[p]):
                # (value, uncertainty)
                lib_props[p] = (row[p], row['u_'+p])
        # isochrones requires fitting a distance to at least one magnitude - use arbitrary value
        lib_props['V'] = 1.

        # create the model and perform the fit
        model = StarModel(dar, **lib_props)
        model.fit(overwrite=True, verbose=False)

        # fill out unknown parameters
        for p in STAR_PROPS:
            value = model.samples[p].quantile(0.5)
            upper_bound = model.samples[p].quantile(MAX_PERCENTILE)
            lower_bound = model.samples[p].quantile(MIN_PERCENTILE)
            if p in lib_props:
                print(p)
                print("\tLibrary values: {0:.2f} +/- {1:.2f}".format(row[p], row['u_'+p]))
                print("\tModel values: {0:.2f}, 2-sigma = ({1:.2f}, {2:.2f})\n".format(
                    value, lower_bound, upper_bound))
                # check for model consistency with known library values
                if (row[p]+row['u_'+p]) < lower_bound or (row[p]-row['u_'+p]) > upper_bound:
                    print("Warning: Model for star {0} had inconsistent values in {1}:".format(
                        row['cps_name'], p))
                    # Save error messages to file
                    if diagnostic:
                        outpath_err = os.path.join(outdir, 'model_errors.txt')
                        with open(outpath_err, 'a') as f:
                            f.write("Inconsistent {0} for star {1}\n".format(p, row['cps_name']))
                            f.write("\tLibrary values: {0:.2f} +/- {1:.2f}\n".format(row[p], row['u_'+p]))
                            f.write("\tModel values: {0:.2f}, 1-sigma = ({1:.2f}, {2:.2f})\n".format(
                            value, lower_bound, upper_bound))

            else:
                # insert unknown values
                stars.loc[i, p] = np.around(value, 2)
                stars.loc[i, 'u_'+p] = np.around(max(upper_bound-value, value-lower_bound), 2)

        # save model
        if diagnostic:
            outpath = os.path.join(outdir, '{0}_model.h5'.format(row['cps_name']))
            model.save_hdf(outpath)

    return stars


def main(catalogdir, cpsdir, outdir, iso, diagnostic, append):
    ### Read in the stars with known stellar parameters and check for those with CPS spectra
    # print("Reading catalogs...")
    # stars, stars_nospectra = read_catalogs(catalogdir, cpsdir)
    # mask = pd.DataFrame()

    # # convert numeric columns
    # for col in STAR_PROPS:
    #     stars[col] = pd.to_numeric(stars[col], errors='coerce')
    #     stars['u_'+col] = pd.to_numeric(stars['u_'+col], errors='coerce')
    #     # Create mask to indicate parameters obtained directly from literature
    #     mask[col] = np.isfinite(stars[col])
    #     mask['u_'+col] = np.isfinite(stars['u_'+col])

    # stars.to_csv(os.path.join(outdir, "libstars_c.csv"))
    # stars_nospectra.to_csv(os.path.join(outdir, "nospectra_c.csv"))
    # mask.to_csv(os.path.join(outdir, "libstars_mask.csv"))

    stars = pd.read_csv(os.path.join(outdir, "libstars_c.csv"), index_col=0)

    ### Fill in remaining parameters with isochrone models
    if iso:
        print("Obtaining isochrone parameters")
        stars = get_isochrone_params(stars, diagnostic=diagnostic, outdir=outdir)

    stars.to_csv(os.path.join(outdir, "libstars_c.csv"))


if __name__ == '__main__':
    psr = ArgumentParser(description="Build the SpecMatch-Emp library from the various catalogs")
    psr.add_argument('catalogdir', type=str, help="Path to catalogs")
    psr.add_argument('cpsdir', type=str, help="Path to CPS spectrum database")
    psr.add_argument('outdir', type=str, help="Path to output directory")
    psr.add_argument('-i', '--isochrones', action='store_true', help="Use isochrones to compute additional parameters")
    psr.add_argument('-d', '--diagnostic', action='store_true', help="Output all intermediate data for diagnostics")
    psr.add_argument('-a', '--append', action='store_true', help="Append to existing library in outdir")
    args = psr.parse_args()

    main(args.catalogdir, args.cpsdir, args.outdir, args.isochrones, args.diagnostic, args.append)
