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

import os
from argparse import ArgumentParser
import re
import warnings

import numpy as np
import pandas as pd

from astroquery.simbad import Simbad
from astropy.io import ascii
from isochrones.dartmouth import Dartmouth_Isochrone
from isochrones import StarModel

from specmatchemp import library
from specmatchemp import shift_spectra
from specmatchemp.io import specmatchio

import time

# relative catalog locations
MANN_FILENAME = "Mann2015/stars.dat"
MANN_README = "Mann2015/ReadMe"
HUBER_FILENAME = "Huber2013/table2.dat"
HUBER_README = "Huber2013/ReadMe"
VONBRAUN_FILENAME = "VonBraun-Figure6.txt"
# BREWER_FILENAME = "spocsii_stars.csv"
BREWER_FILENAME = "brewer_cut.csv"
CPS_INDEX = "cps_templates.csv"
CPS_SPECTRA_DIR = "iodfitsdb/"

STAR_PROPS = ['Teff', 'radius', 'logg', 'feh', 'mass', 'age']
NOSPECTRA_COLS = ['name', 'source']

MIN_PERCENTILE = 0.16
MAX_PERCENTILE = 0.84

WAVLIM = (5000, 6400)

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
    elif starname.upper().startswith(('BD+', 'EPIC', 'HTR', 'HATS', 'HIP', 'HII')):
        starname = starname.strip('*')
        cps_search_str = starname
    # don't match any other identifiers
    else:
        cps_search_str = 'NOSTARFOUND'

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
    print("\tReading Brewer catalog")
    brewer_stars, brewer_nospec = read_brewer(catalogdir, cps_list)
    stars = pd.concat((stars, brewer_stars), ignore_index=True)
    stars_nospectra = stars_nospectra.append(brewer_nospec)

    print("\tReading Mann catalog")
    mann_stars, mann_nospec = read_mann(catalogdir, cps_list)
    stars = pd.concat((stars, mann_stars), ignore_index=True)
    stars_nospectra = stars_nospectra.append(mann_nospec)

    print("\tReading von Braun catalog")
    vonbraun_stars, vonbraun_nospec = read_vonbraun(catalogdir, cps_list)
    stars = pd.concat((stars, vonbraun_stars), ignore_index=True)
    stars_nospectra = stars_nospectra.append(vonbraun_nospec)

    print("\tReading Huber catalog")
    huber_stars, huber_nospec = read_huber(catalogdir, cps_list)
    stars = pd.concat((stars, huber_stars), ignore_index=True)
    stars_nospectra = stars_nospectra.append(huber_nospec)

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
        model = StarModel(dar, use_emcee=True, **lib_props)
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
                    print("\tModel values: {0:.2f}, 1-sigma = ({1:.2f}, {2:.2f})\n".format(
                        value, lower_bound, upper_bound))
                    # Save error messages to file
                    if diagnostic:
                        outpath_err = os.path.join(outdir, 'model_errors.txt')
                        with open(outpath_err, 'a') as f:
                            f.write("Inconsistent {0} for star {1}".format(p, row['cps_name']))
                            f.write("\tLibrary values: {0:.2f} +/- {1:.2f}".format(row[p], row['u_'+p]))
                            f.write("\tModel values: {0:.2f}, 1-sigma = ({1:.2f}, {2:.2f})".format(
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

def shift_library(stars, cpsdir, shift_reference, diagnostic=False, outdir='~/'):
    """Reads in spectra and shifts them to a reference spectrum and interpolates onto
    a constant log-lambda scale.

    Args:
        stars (pd.DataFrame) : star library
        shift_reference (pd.DataFrame) : DataFrame containing a list of reference spectra.
            This allows bootstrapping in order to shift spectra in different regions of the
            parameter space. 
            The DataFrame should have the following columns:
                index: Determines order of spectra to shift
                min_t, max_t: Teff range for this reference spectrum.
                ref_spectrum: A path to a FITS file, or the obs value of a spectrum that
                    is in the library and has previously been shifted.
        diagnostic (bool)   : whether to save the diagnostic data
        outdir (str)        : directory to save diagnostic data

    Returns:
        stars (pd.DataFrame) : Updated star library with filled lib_obs, lib_index columns
        wav (np.ndarray)     : 1-dimensional array containing the common wavelength scale.
        spectra(np.ndarray)  : 3-dimensional array containing the spectra, u_spectra,
            indexed by the star library index.
    """
    # obey the specified shifting order
    wav = None
    spectra = None
    shift_reference = shift_reference.sort_index()

    stars.is_copy=False

    num_spec = len(stars)
    current_spec = 1

    for i in range(len(shift_reference)):
        min_t = shift_reference.iloc[i].min_t
        max_t = shift_reference.iloc[i].max_t
        ref_spectrum = shift_reference.iloc[i].ref_spectrum

        # for the first round, we must specify a file (or previous library spectrum)
        if i == 0:
            if not os.path.isfile(ref_spectrum):
                print("Error: Unable to find {0} as reference.".format(ref_spectrum))
                return stars, None, None
            w_ref, s_ref, serr_ref, hdr_ref = specmatchio.read_standard_spectrum(ref_spectrum, (4900,6500))
            # use this reference wavelength scale for the library
            wav, s, serr = specmatchio.truncate_spectrum(WAVLIM, w_ref, s_ref, serr_ref)
            spectra = np.empty((0,2,len(wav)))
        # in other rounds, we can look for a file or a previously-shifted spectrum for bootstrapping
        else:
            # search for file
            if os.path.isfile(ref_spectrum):
                w_ref, s_ref, serr_ref, hdr_ref = specmatchio.read_standard_spectrum(ref_spectrum, (4900,6500))
            # search in previously-shifted spectra
            else:
                pattern = '^' + ref_spectrum + '$'
                # search previously-shifted stars
                if not stars[stars.lib_obs.str.match(pattern)].empty:
                    ref_idx = int(stars[stars.lib_obs.str.match(pattern)].iloc[0].lib_index)
                    s_ref = spectra[ref_idx,0]
                    serr_ref = spectra[ref_idx, 1]
                    w_ref = wav
                else:
                    print("Error: Unable to find {0} as reference.".format(ref_spectrum))
                    continue

        # now shift spectra in each group
        query = '{0:.0f} <= Teff < {1:.0f}'.format(min_t, max_t)
        stars_grouped = stars.query(query)
        for targ_idx, targ_params in stars_grouped.iterrows():
            print("Shifting star {0:d} of {1:d}".format(current_spec, num_spec))
            current_spec+=1
            # find spectrum for the star
            obs_list = re.sub("[\[\]\'']", '', targ_params.obs).split()
            specpath = None
            ## loop through list of observations
            for obs in obs_list:
                specpath = os.path.join(cpsdir, CPS_SPECTRA_DIR, obs+'.fits')
                if os.path.isfile(specpath):
                    # record which observation was used
                    stars.loc[targ_idx, 'lib_obs'] = obs
                    break
                else:
                    specpath=None
            if specpath is None:
                print("Could not find any spectra for star {0}".format(targ_params.cps_name))
                continue

            try:
                ## read in spectrum
                w_targ, s_targ, serr_targ, hdr_targ = specmatchio.read_hires_spectrum(specpath)

                # shift spectrum
                outfile = os.path.join(outdir,'/shift_data/{0}.txt'.format(stars.loc[targ_idx].lib_obs))
                diag_hdr = '# Star: {0}\n# Reference: {1}\n'.format(targ_params.cps_name, ref_spectrum)

                s_adj, serr_adj, w_adj = shift_spectra.adjust_spectra(s_targ, serr_targ, w_targ,\
                    s_ref, serr_ref, w_ref, diagnostic=diagnostic, outfile=outfile, diagnostic_hdr=diag_hdr)

                # flatten spectrum within limits
                w_flat, s_flat, serr_flat = shift_spectra.flatten(w_adj, s_adj, serr_adj, w_ref=wav, wavlim=WAVLIM)

            except Exception as e:
                print("Error: Failed to shift star {0}, spectrum {1}".format(targ_params.cps_name, targ_params.lib_obs))
                print("{0}".format(e))
                continue


            # append spectrum to spectra table
            stars.loc[targ_idx, 'lib_index'] = len(spectra)
            spectra = np.vstack((spectra, [[s_flat, serr_flat]]))

            # calculate the signal-to-noise-ratio
            stars.loc[targ_idx, 'snr'] = np.nanpercentile(1/serr_flat, 90)

    # eliminate stars with no spectra
    stars = stars[np.logical_not(np.isnan(stars.lib_index))]
    return stars, wav, spectra

def main(catalogdir, cpsdir, shift_reference_path, outdir, diagnostic, append):
    ### 1. Read in the stars with known stellar parameters and check for those with CPS spectra
    # print("Step 1: Reading catalogs...")
    # stars, stars_nospectra = read_catalogs(catalogdir, cpsdir)
    # # convert numeric columns
    # for col in STAR_PROPS:
    #     stars[col] = pd.to_numeric(stars[col], errors='coerce')
    #     stars['u_'+col] = pd.to_numeric(stars['u_'+col], errors='coerce')

    # ### 2. Use isochrones package to obtain the remaining, unknown stellar parameters
    # print("Step 2: Obtaining isochrone parameters...")
    # stars = get_isochrone_params(stars, diagnostic=diagnostic, outdir=outdir)

    # stars.to_csv(os.path.join(outdir, "libstars.csv"))
    # stars_nospectra.to_csv(os.path.join(outdir, "stars_nospectra.csv"))

    # stars['obs'] = stars['obs'].astype(str)

    ################################################################
    stars = pd.read_csv("./lib/libstars.csv", index_col=0)
    stars['lib_obs'] = stars['lib_obs'].astype(str)
    ################################################################

    stars.reset_index(drop=True,inplace=True)
    stars = stars.loc[[234,]]

    ### 3. Shift library spectra onto a constant log-lambda scale
    print("Step 3: Shifting library spectra...")
    shift_ref = pd.read_csv(shift_reference_path, index_col=0)
    stars, wav, spectra = shift_library(stars, cpsdir, shift_ref, diagnostic=diagnostic, outdir=outdir)

    stars.to_csv(os.path.join(outdir, "libstars_shifted.csv"))

    ### 4. Create and save the library
    print("Step 4: Saving library...")
    stars = stars.drop('obs', axis=1)
    lib = library.Library(wav, spectra, stars, wavlim=WAVLIM)
    lib.to_hdf('./lib/library.h5')


if __name__ == '__main__':
    psr = ArgumentParser(description="Build the SpecMatch-Emp library from the various catalogs")
    psr.add_argument('catalogdir', type=str, help="Path to catalogs")
    psr.add_argument('cpsdir', type=str, help="Path to CPS spectrum database")
    psr.add_argument('shift_reference', type=str, help="Path to spectrum shifting reference list")
    psr.add_argument('outdir', type=str, help="Path to output directory")
    psr.add_argument('-d', '--diagnostic', action='store_true', help="Output all intermediate data for diagnostics")
    psr.add_argument('-a', '--append', action='store_true', help="Append to existing library in outdir")
    args = psr.parse_args()

    main(args.catalogdir, args.cpsdir, args.shift_reference, args.outdir, args.diagnostic, args.append)
