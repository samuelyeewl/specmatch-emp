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
import time

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import multiprocessing as mp

from astroquery.simbad import Simbad
from astropy.io import ascii
from isochrones.dartmouth import Dartmouth_Isochrone
from isochrones import StarModel

from specmatchemp import library

# relative catalog locations
MANN_FILENAME = "Mann2015/stars.dat"
MANN_README = "Mann2015/ReadMe"
HUBER_FILENAME = "Huber2013/table2.dat"
HUBER_README = "Huber2013/ReadMe"
VONBRAUN_FILENAME = "VonBraun-Figure6.txt"
BREWER_FILENAME = "spocsii_stars.csv"
CPS_INDEX = "cps_templates.csv"
CPS_SPECTRA_DIR = "iofitsdb/"

STAR_PROPS = ['Teff', 'radius', 'logg', 'feh', 'mass', 'age']
NOSPECTRA_COLS = ['name', 'source']

MIN_PERCENTILE = 0.16
MAX_PERCENTILE = 0.84

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
        num = int(starname.split('KOI')[1].strip('-'))
        cps_search_str = 'K'+'{0:05d}'.format(num)+'|CK'+'{0:05d}'.format(num)
        
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

    cps_search_str = '^' + cps_search_str + '$'
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

    stars = pd.DataFrame(columns=library.LIB_COLS)
    nospectra = pd.DataFrame(columns=NOSPECTRA_COLS)

    for row in brewer_data:
        try:
            query_result = find_spectra(row['NAME'], cps_list)
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

    stars = pd.DataFrame(columns=library.LIB_COLS)
    nospectra = pd.DataFrame(columns=NOSPECTRA_COLS)

    for row in mann_data:
        try:
            query_result = find_spectra(row['CNS3'], cps_list)
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

    stars = pd.DataFrame(columns=library.LIB_COLS)
    nospectra = pd.DataFrame(columns=NOSPECTRA_COLS)

    for row in vb_data:
        try:
            query_result = find_spectra(row['Star'], cps_list)
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

    stars = pd.DataFrame(columns=library.LIB_COLS)
    nospectra = pd.DataFrame(columns=NOSPECTRA_COLS)

    for row in huber_data:
        try:
            query_result = find_spectra('KOI'+str(row['KOI']), cps_list)
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
    stars = pd.DataFrame(columns=library.LIB_COLS)

    # Create dataframe to store stars which have no spectra
    stars_nospectra = pd.DataFrame(columns=NOSPECTRA_COLS)

    # Read catalogs
    # brewer_stars, brewer_nospec = read_brewer(catalogdir, cps_list)
    # stars = stars.append(brewer_stars)
    # stars_nospectra = stars_nospectra.append(brewer_nospec)

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

def get_isochrone_params(stars, diagnostic=False, outdir='~/'):
    """Fill out parameter table with values obtained from isochrone package

    Args:
        stars (pd.DataFrame): star library
        diagnostic (bool)   : whether to save the fitted model as a file
        outdir (str)        : directory to save fitted models
    Returns:
        stars (pd.DataFrame): star library with updated parameters
    """
    dar = Dartmouth_Isochrone()

    for i, row in stars.iterrows():
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
        model.fit_mcmc()

        # fill out unknown parameters
        for p in STAR_PROPS:
            value = model.samples[p].quantile(0.5)
            upper_bound = model.samples[p].quantile(MAX_PERCENTILE)
            lower_bound = model.samples[p].quantile(MIN_PERCENTILE)
            if p in lib_props:
                # check for model consistency with known library values
                if (row[p]+row['u_'+p]) < lower_bound or (row[p]-row['u_'+p]) > upper_bound:
                    print("Warning: Model for star {0} had inconsistent values in {1}:".format(
                        row['cps_name'], p))
                    print("\tLibrary values: {0:.2f} +/- {1:.2f}".format(row[p], row['u_'+p]))
                    print("\tModel values: {0:.2f}, 1-sigma = ({1:.2f}, {2:.2f})".format(
                        value, lower_bound, upper_bound))
            else:
                # insert unknown values
                stars.loc[i, p] = np.around(value, 2)
                stars.loc[i, 'u_'+p] = np.around(max(upper_bound-value, value-lower_bound), 2)

        # save model
        if diagnostic:
            outpath = os.path.join(outdir, 'isochrone_models/')
            if not os.path.exists(outpath):
                os.makedirs(outpath)
            outpath = os.path.join(outpath, '{0}_model.h5'.format(row['cps_name']))
            model.save_hdf(outpath)

    return stars

def shift_library(stars, diagnostic=False, outdir='~/'):

    return

def main(catalogdir, cpsdir, outdir, diagnostic):
    ### 1. Read in the stars with known stellar parameters and check for those with CPS spectra
    # stars, stars_nospectra = read_catalogs(catalogdir, cpsdir)
    # # convert numeric columns
    # for col in STAR_PROPS:
    #     stars[col] = pd.to_numeric(stars[col], errors='coerce')
    #     stars['u_'+col] = pd.to_numeric(stars['u_'+col], errors='coerce')
    # stars.to_csv(os.path.join(outdir, "libstars_small.csv"))
    # stars_nospectra.to_csv(os.path.join(outdir, "stars_nospectra.csv"))
    ################################################################
    stars = pd.read_csv("./lib/libstars_small.csv", index_col=0)
    ################################################################

    stars = stars.head(8)

    start = time.time()

    ### 2. Use isochrones package to obtain the remaining, unknown stellar parameters
    stars = get_isochrone_params(stars, diagnostic=diagnostic, outdir=outdir)

    end = time.time()

    print("time to fit 8 stars = {0:d}".format(end-start))

    ### 3. Shift library spectra onto a constant log-lambda scale
    # stars = shift_library(stars, diagnostic=diagnostic, outdir=outdir)





if __name__ == '__main__':
    psr = ArgumentParser(description="Build the SpecMatch-Emp library from the various catalogs")
    psr.add_argument('catalogdir', type=str, help="Path to catalogs")
    psr.add_argument('cpsdir', type=str, help="Path to CPS spectrum database")
    psr.add_argument('outdir', type=str, help="Path to output directory")
    psr.add_argument('-d', '--diagnostic', action='store_true', help="Output all intermediate data for diagnostics")
    args = psr.parse_args()

    main(args.catalogdir, args.cpsdir, args.outdir, args.diagnostic)
