#!/usr/bin/env python
"""
@filename build_library.py

Builds the SpecMatch-Emp library from the Huber, Mann, Von Braun and Brewer catalogs.
1.  Reads in the stars with known stellar parameters and checks for spectra in the CPS
    spectrum database.
2.  Uses isochrones package to obtain a full set of stellar parameters (radius, mass, logg)
    from those which are known.
3.  Shifts library spectra onto a constant log-lambda scale.
4.  Saves library table and spectra as a Library object in a hdf file.
"""

import glob
import os
from argparse import ArgumentParser

import re
import warnings

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astroquery.simbad import Simbad
from astropy.io import ascii

# relative catalog locations
MANN_FILENAME = "Mann2015/stars.dat"
MANN_README = "Mann2015/ReadMe"
HUBER_FILENAME = "Huber2013/table2.dat"
HUBER_README = "Huber2013/ReadMe"
VONBRAUN_FILENAME = "VonBraun-Figure6.txt"
BREWER_FILENAME = "spocsii_stars.csv"
CPS_INDEX = "cps_templates.csv"
CPS_SPECTRA_DIR = "iofitsdb/"

LIB_COLS = ['cps_name', 'obs', 'Teff', 'u_Teff', 'radius', 'u_radius', 'logg', 'u_logg','FeH', 'u_FeH',
           'mass', 'u_mass', 'vsini', 'source', 'source_name']
NOSPECTRA_COLS = ['name', 'source']

def check_cps_database(starname, cps_list):
    """Checks the CPS databse if a spectrum is available for the given identifier.

    Parses the provided identifier and converts it to the format used in the CPS library.

    Args:
        starname (str): an identifier for the star
        cps_list (pd.DataFrame): Pandas dataframe with list of stars with spectra
    
    Returns:
        Subset of cps_list corresponding to the observations for the given star.
    """
    starname = re.sub('\s', '', starname)
    cps_search_str = ''
    
    # HD identifiers appear as just the number in the CPS database
    if starname.startswith('HD'):
        cps_search_str = '^'+starname.split('HD')[1]+'$'
    
    # Gliese database identifiers may appear as GJ or GL
    elif starname.startswith('GJ'):
        num = starname.split('GJ')[1]
        cps_search_str = 'GJ'+str(num)+'|GL'+str(num)
    elif starname.startswith('Gj'):
        num = starname.split('Gj')[1]
        cps_search_str = 'GJ'+str(num)+'|GL'+str(num)
    elif starname.startswith('GL'):
        num = starname.split('GL')[1]
        cps_search_str = 'GJ'+str(num)+'|GL'+str(num)
    elif starname.startswith('Gl'):
        num = starname.split('Gl')[1]
        cps_search_str = 'GJ'+str(num)+'|GL'+str(num)
    
    # KIC identifiers should not have dashes and may be have leading zeroes
    elif starname.startswith('KIC'):
        num = re.sub('[^0-9]', '', starname).lstrip('0')
        cps_search_str = 'KIC'+'0*'+str(num)
    
    # KOI identifiers may be prefixed by K or CK and have 5 digits padded by leading zeros
    elif starname.startswith('KOI'):
        num = starname.split('KOI')[1].strip('-')
        cps_search_str = 'K'+'0*'+str(num)+'|CK'+'0*'+str(num)
        
    # WASP identifiers may be either WASP-x or WASPx, x is the idnum
    elif starname.startswith('WASP'):
        num = re.sub('[^0-9]', '', starname)
        cps_search_str = 'WASP'+'[\-]*'+str(num)

    # COROT identifiers may be either COROT-x or COROTx
    elif starname.upper().startswith('COROT'):
        num = re.sub('[^0-9]', '', starname)
        cps_search_str = 'COROT'+'[\-]*'+str(num)
    
    # TRES identifiers may be either TRES-x or TRESx
    elif starname.upper().startswith('TRES'):
        num = re.sub('[^0-9]', '', starname)
        cps_search_str = 'TRES'+'[\-]*'+str(num)
    
    # else the search string is simply the name
    # e.g. BD+, EPIC, HTR, HATS, HIP, HII
    else:
        starname = starname.strip('*')
        cps_search_str = starname

    return cps_list[cps_list['name'].str.match(cps_search_str)]

def find_spectra(starname, cps_list):
    """Checks the CPS database for the star with the provided identifier.
    If not found, queries SIMBAD database for alternative idenitifiers to
    check.

    Args:
        starname (str): an identifier for the star
        cps_list (pd.DataFrame): Pandas dataframe with list of stars with spectra
    
    Returns:
        Subset of cps_list corresponding to the observations for the given star.
    """

    result = check_cps_database(starname, cps_list)
    
    if result.empty:
        # ignore warnings from Simbad query - we check if there are any results
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            s_query_result = Simbad.query_objectids(starname)
            if s_query_result is not None:
                for r in s_query_result:
                    result = check_cps_database(r['ID'].decode('ascii'), cps_list)
                    if not result.empty:
                        break
                    
    return result

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

    for row in brewer_data:
        try:
            query_result = find_spectra(row['NAME'])
            if not query_result.empty:
                new_row = {}
                new_row['cps_name'] = query_result.iloc[0].name
                new_row['obs'] = query_result.obs.values
                new_row['Teff'] = row['TEFF']
                new_row['u_Teff'] = 25
                new_row['logg'] = row['LOGG']
                new_row['u_logg'] = 0.028
                new_row['FeH'] = row['FEH']
                new_row['u_FeH'] = 0.010
                new_row['vsini'] = row['VSINI']
                new_row['source'] = 'Brewer'
                new_row['source_name'] = row['NAME']
                stars = stars.append(pd.Series(new_row), ignore_index=True)
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
            query_result = find_spectra(row['CNS3'])
            if not query_result.empty:
                new_row = {}
                new_row['cps_name'] = query_result.iloc[0].name
                new_row['obs'] = query_result.obs.values
                new_row['Teff'] = row['Teff']
                new_row['u_Teff'] = row['e_Teff']
                new_row['radius'] = row['R']
                new_row['u_radius'] = row['e_R']
                new_row['FeH'] = row['[Fe/H]']
                new_row['u_FeH'] = row['e_[Fe/H]']
                new_row['mass'] = row['M']
                new_row['u_mass'] = row['e_M']
                new_row['source'] = 'Mann'
                new_row['source_name'] = row['CNS3']
                stars = stars.append(pd.Series(new_row), ignore_index=True)
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
            query_result = find_spectra(row['Star'])
            if not query_result.empty:
                new_row = {}
                new_row['cps_name'] = query_result.iloc[0].name
                new_row['obs'] = query_result.obs.values
                new_row['Teff'] = row['Teff']
                new_row['u_Teff'] = row['eTeff']
                new_row['radius'] = row['Radius']
                new_row['u_radius'] = row['eRadius']
                new_row['FeH'] = row['FeH']
                new_row['source'] = 'Von Braun'
                new_row['source_name'] = row['Star']
                stars = stars.append(pd.Series(new_row), ignore_index=True)
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

    stars = pd.DataFrame(columns=LIB_COLS)
    nospectra = pd.DataFrame(columns=NOSPECTRA_COLS)

    for row in huber_data:
        try:
            query_result = find_spectra('KOI'+str(row['KOI']), cps_list)
            if not query_result.empty:
                new_row = {}
                new_row['cps_name'] = query_result.iloc[0].name
                new_row['obs'] = query_result.obs.values
                new_row['Teff'] = row['Teff']
                new_row['u_Teff'] = row['e_Teff']
                new_row['radius'] = row['Rad']
                new_row['u_radius'] = row['e_Rad']
                new_row['FeH'] = row['[Fe/H]']
                new_row['u_FeH'] = row['e_[Fe/H]']
                new_row['mass'] = row['Mass']
                new_row['u_mass'] = row['e_Mass']
                new_row['source'] = 'Huber'
                new_row['source_name'] = 'KOI'+str(row['KOI'])
                stars = stars.append(pd.Series(new_row), ignore_index=True)
        except:
            new_row = {}
            new_row['name'] = 'KOI'+str(row['KOI'])
            new_row['source'] = 'Huber'
            nospectra = nospectra.append(pd.Series(new_row), ignore_index=True)

    return stars, nospectra

def read_catalogs(catalogdir, cpsdir):
    """Reads in the Brewer, Mann, Von Braun and Huber catalogs

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
    brewer_stars, brewer_nospec = read_brewer(catalogdir, cps_list)
    stars = stars.append(brewer_stars)
    stars_nospectra = stars_nospectra.append(brewer_nospec)

    mann_stars, mann_nospec = read_mann(catalogdir, cps_list)
    stars = stars.append(mann_stars)
    stars_nospectra = stars_nospectra.append(mann_nospec)

    vonbraun_stars, vonbraun_nospec = read_vonbraun(catalogdir, cps_list)
    stars = stars.append(vonbraun_stars)
    stars_nospectra = stars_nospectra.append(vonbraun_nospec)

    huber_stars, huber_nospec = read_huber(catalogdir, cps_list)
    stars = stars.append(huber_stars)
    stars_nospectra = stars_nospectra.append(huber_nospec)

    return stars, stars_nospectra

def main(catalogdir, cpsdir, outdir):
    ### 1. Read in the stars with known stellar parameters and check for those with CPS spectra
    stars, stars_nospectra = read_catalogs(catalogdir, cpsdir)
    stars.to_csv(os.path.join(outdir, "libstars.csv"))
    stars_nospectra.to_csv(os.path.join(outdir, "stars_nospectra.csv"))



if __name__ == '__main__':
    psr = ArgumentParser(description="Build the SpecMatch-Emp library from the various catalogs")
    psr.add_argument('catalogdir', type=str, help="Path to catalogs")
    psr.add_argument('cpsdir', type=str, help="Path to CPS spectrum database")
    psr.add_argument('outdir', type=str, help="Path to output directory")
    args = psr.parse_args()

    main(args.catalogdir, args.cpsdir, args.outdir)
