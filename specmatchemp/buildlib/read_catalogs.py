#!/usr/bin/env python
"""
@filename read_catalogs.py

Builds the parameter table of the SpecMatch-Emp library by reading in the
various catalogs, checking for spectra in the CPS database.
Saves the parameters as a DataFrame in a csv file.
"""

from __future__ import print_function

import os
from argparse import ArgumentParser

import pandas as pd

from astropy.io import ascii

from specmatchemp.library import Library
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
KDWARF_FILENAME = "kdwarfs_cut.csv"
CPS_INDEX = "cps_templates.csv"
CPS_SPECTRA_DIR = "iodfitsdb/"

NOSPECTRA_COLS = ['name', 'source']


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

    stars = pd.DataFrame(columns=Library.LIB_COLS)
    nospectra = pd.DataFrame(columns=NOSPECTRA_COLS)

    for row in brewer_data:
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
    mann_data = ascii.read(os.path.join(catalogdir, MANN_FILENAME),
                           readme=os.path.join(catalogdir, MANN_README))

    mann_data = mann_data[~mann_data['CNS3'].mask]

    stars = pd.DataFrame(columns=Library.LIB_COLS)
    nospectra = pd.DataFrame(columns=NOSPECTRA_COLS)

    for row in mann_data:
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
            logg, u_logg = utils.calc_logg(row['R'], row['e_R'],
                                           row['M'], row['e_M'])
            new_row['logg'] = logg
            new_row['u_logg'] = u_logg

            stars = stars.append(pd.Series(new_row), ignore_index=True)
        else:
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

    stars = pd.DataFrame(columns=Library.LIB_COLS)
    nospectra = pd.DataFrame(columns=NOSPECTRA_COLS)

    for row in vb_data:
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
    huber_data = ascii.read(os.path.join(catalogdir, HUBER_FILENAME),
                            readme=os.path.join(catalogdir, HUBER_README))
    huber_data = huber_data[huber_data['f_KOI'] != '*']

    stars = pd.DataFrame(columns=Library.LIB_COLS)
    nospectra = pd.DataFrame(columns=NOSPECTRA_COLS)

    for row in huber_data:
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
            logg, u_logg = utils.calc_logg(row['Rad'], row['e_Rad'],
                                           row['Mass'], row['e_Mass'])
            new_row['logg'] = logg
            new_row['u_logg'] = u_logg

            stars = stars.append(pd.Series(new_row), ignore_index=True)
        else:
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

    stars = pd.DataFrame(columns=Library.LIB_COLS)
    nospectra = pd.DataFrame(columns=NOSPECTRA_COLS)

    for idx, row in ramirez_data.iterrows():
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

    stars = pd.DataFrame(columns=Library.LIB_COLS)
    nospectra = pd.DataFrame(columns=NOSPECTRA_COLS)

    for idx, row in c_data.iterrows():
        query_result = cpsutils.find_spectra(row['Name'], cps_list)
        if not query_result.empty:
            # Calculate stellar radius from angular diameter and parallax
            radius, u_radius = utils.calc_radius(row['Plx'], row['e_Plx'],
                                                 row['Diam'], row['e_Diam'])
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
    b_data = ascii.read(os.path.join(catalogdir, BRUNTT_FILENAME),
                        readme=os.path.join(catalogdir, BRUNTT_README))

    stars = pd.DataFrame(columns=Library.LIB_COLS)
    nospectra = pd.DataFrame(columns=NOSPECTRA_COLS)

    for row in b_data:
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

    return stars, nospectra


def read_kdwarfs(catalogdir, cps_list):
    """Read in K-Dwarf catalog"""
    kdwarfs = pd.read_csv(os.path.join(catalogdir, KDWARF_FILENAME))

    stars = pd.DataFrame(columns=Library.LIB_COLS)
    nospectra = pd.DataFrame(columns=NOSPECTRA_COLS)

    for idx, row in kdwarfs.iterrows():
        query_result = cpsutils.find_spectra(row['name'], cps_list)
        if not query_result.empty:
            new_row = {}
            new_row['cps_name'] = str(query_result.iloc[0]['name'])
            new_row['obs'] = query_result.obs.values
            new_row['Teff'] = row['teff_derived']
            new_row['u_Teff'] = row['e_teff_derived']
            new_row['feh'] = row['fe']
            new_row['u_feh'] = 0.1
            new_row['logg'] = row['logg']
            new_row['u_logg'] = row['u_logg']
            new_row['mass'] = row['mass']
            new_row['u_mass'] = row['u_mass']
            new_row['radius'] = row['radius']
            new_row['u_radius'] = row['u_radius']
            new_row['age'] = row['age']
            new_row['u_age'] = row['u_age']
            new_row['source'] = 'Gaidos'
            new_row['source_name'] = row['name']

            stars = stars.append(pd.Series(new_row), ignore_index=True)
        else:
            new_row = {}
            new_row['name'] = row['name']
            new_row['source'] = 'Gaidos'
            nospectra = nospectra.append(pd.Series(new_row), ignore_index=True)

    return stars, nospectra


def read_catalogs(catalogdir, cpspath):
    """Reads in the catalogs

    Args:
        catalogdir (str): Path to catalog directory
        cpsdir (str): Path to CPS list

    Returns:
        stars (pd.DataFrame): Stars in source which have CPS spectra
        nospec (pd.DataFrame): Stars in source which don't have CPS spectra
    """
    # Read in list of CPS spectra
    cps_list = pd.read_csv(cpspath)

    # Create dataframe to store found stars
    stars = pd.DataFrame(columns=Library.LIB_COLS)

    # Create dataframe to store stars which have no spectra
    stars_nospectra = pd.DataFrame(columns=NOSPECTRA_COLS)

    # Read catalogs
    print("Reading Brewer catalog")
    brewer_stars, brewer_nospec = read_brewer(catalogdir, cps_list)
    stars = pd.concat((stars, brewer_stars), ignore_index=True)
    stars_nospectra = stars_nospectra.append(brewer_nospec)
    print("\t{0:d} stars with spectra from Brewer catalog"
          .format(len(brewer_stars)))

    print("Reading Mann catalog")
    mann_stars, mann_nospec = read_mann(catalogdir, cps_list)
    stars = pd.concat((stars, mann_stars), ignore_index=True)
    stars_nospectra = stars_nospectra.append(mann_nospec)
    print("\t{0:d} stars with spectra from Mann catalog"
          .format(len(mann_stars)))

    print("Reading von Braun catalog")
    vonbraun_stars, vonbraun_nospec = read_vonbraun(catalogdir, cps_list)
    stars = pd.concat((stars, vonbraun_stars), ignore_index=True)
    stars_nospectra = stars_nospectra.append(vonbraun_nospec)
    print("\t{0:d} stars with spectra from von Braun catalog"
          .format(len(vonbraun_stars)))

    # print("Reading Huber catalog")
    # huber_stars, huber_nospec = read_huber(catalogdir, cps_list)
    # stars = pd.concat((stars, huber_stars), ignore_index=True)
    # stars_nospectra = stars_nospectra.append(huber_nospec)
    # print("\t{0:d} stars with spectra from Huber catalog"
    #       .format(len(huber_stars)))

    # print("Reading Ramirez catalog")
    # ramirez_stars, ramirez_nospec = read_ramirez(catalogdir, cps_list)
    # stars = pd.concat((stars, ramirez_stars), ignore_index=True)
    # stars_nospectra = stars_nospectra.append(ramirez_nospec)
    # print("\t{0:d} stars with spectra from Ramirez catalog"
    #       .format(len(ramirez_stars)))

    # print("Reading Casagrande catalog")
    # c_stars, c_nospec = read_casagrande(catalogdir, cps_list)
    # stars = pd.concat((stars, c_stars), ignore_index=True)
    # stars_nospectra = stars_nospectra.append(c_nospec)
    # print("\t{0:d} stars with spectra from Casagrande catalog"
    #       .format(len(c_stars)))

    print("Reading Bruntt catalog")
    b_stars, b_nospec = read_bruntt(catalogdir, cps_list)
    stars = pd.concat((stars, b_stars), ignore_index=True)
    stars_nospectra = stars_nospectra.append(b_nospec)
    print("\t{0:d} stars with spectra from Bruntt catalog"
          .format(len(b_stars)))

    print("Reading K Dwarfs catalog")
    k_stars, k_nospec = read_kdwarfs(catalogdir, cps_list)
    stars = pd.concat((stars, k_stars), ignore_index=True)
    stars_nospectra = stars_nospectra.append(k_nospec)
    print("\t{0:d} stars with spectra from K Dwarf catalog"
          .format(len(b_stars)))

    dups = stars[stars.duplicated(subset='cps_name', keep=False)].sort_values(
        by='cps_name')
    dups_vb = dups[~dups.source.str.contains('Von Braun|Bruntt')]
    idxs = dups_vb.index
    stars.drop(idxs, inplace=True)

    print("Removing {0:d} duplicates, favoring von Braun and Bruntt data"
          .format(len(idxs)))

    idxs = stars.query('Teff > 7000 | feh < -1.0').index
    stars.drop(idxs, inplace=True)

    print("Removing {0:d} stars outside our parameter space"
          .format(len(idxs)))

    print("Total of {0:d} stars read".format(len(stars)))

    return stars, stars_nospectra


def main(catalogdir, cpspath, outdir, append):
    """Read in the stars with known stellar parameters and check for those
    with CPS spectra"""
    print("Reading catalogs...")
    stars, stars_nospectra = read_catalogs(catalogdir, cpspath)

    # convert numeric columns
    for col in Library.STAR_PROPS:
        stars[col] = pd.to_numeric(stars[col], errors='coerce')
        stars['u_'+col] = pd.to_numeric(stars['u_'+col], errors='coerce')

    if append:
        mode = 'a'
        header = False
    else:
        mode = 'w'
        header = True

    stars.reset_index(drop=True, inplace=True)
    stars_nospectra.reset_index(drop=True, inplace=True)
    stars.to_csv(os.path.join(outdir, "libstars.csv"), mode=mode,
                 header=header)
    stars_nospectra.to_csv(os.path.join(outdir, "nospectra.csv"), mode=mode,
                           header=header)


if __name__ == '__main__':
    psr = ArgumentParser(description="Read catalogs to produce library")
    psr.add_argument('catalogdir', type=str, help="Path to catalogs")
    psr.add_argument('cpslist', type=str, help="Path to CPS spectrum list")
    psr.add_argument('outdir', type=str, help="Path to output directory")
    psr.add_argument('-a', '--append', action='store_true',
                     help="Append to existing library in outdir")
    args = psr.parse_args()

    main(args.catalogdir, args.cpslist, args.outdir, args.append)
