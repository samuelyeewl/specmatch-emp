"""
@filename buildlib/cpsutils.py

Helper functions to find spectra in the CPS database
"""
import os
import re
import warnings
import numpy as np

from specmatchemp import SPECMATCHDIR
from astroquery.simbad import Simbad
from astropy.io import fits


def check_cps_database(starname, cps_list):
    """Checks the CPS databse if a spectrum is available for the given
    identifier.

    Parses the provided identifier and converts it to the format used in the
    CPS library.

    Args:
        starname (str): an identifier for the star
        cps_list (pd.DataFrame): Pandas dataframe with list of stars with
            spectra

    Returns:
        Subset of cps_list corresponding to the observations for the given
        star.
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

    # KOI identifiers may be prefixed by K or CK and have 5 digits padded by
    # leading zeros
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
    elif starname.upper().startswith(('BD+', 'EPIC', 'HTR',
                                      'HATS', 'HIP', 'HII')):
        # starname = starname.strip('*')
        # cps_search_str = starname
        cps_search_str = 'NOSTARFOUND'
    # don't match any other identifiers
    else:
        cps_search_str = 'NOSTARFOUND'

    cps_search_str = '^' + cps_search_str + '$'
    return cps_list[cps_list['name'].str.match(cps_search_str)].copy()


def find_spectra(starname, cps_list, specdir=""):
    """Checks the CPS database for the star with the provided identifier.
    If not found, queries SIMBAD database for alternative idenitifiers to
    check.

    Args:
        starname (str): an identifier for the star
        cps_list (pd.DataFrame): Pandas dataframe with list of stars with
            spectra

    Returns:
        pd.DataFrame: Subset of cps_list corresponding to the observations for
        the given star.
    """
    specdir = "/Users/samuel/Dropbox/SpecMatch-Emp/spectra/iodfitsdb/"
    result = check_cps_database(starname, cps_list)

    if result.empty:
        # ignore warnings from Simbad query - we check if there are any results
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            s_query_result = Simbad.query_objectids(starname)
            if s_query_result is not None:
                for r in s_query_result:
                    result = check_cps_database(r['ID'].decode('ascii'),
                                                cps_list)
                    if not result.empty:
                        break

    # Check SNR for each spectrum and rank results
    if not result.empty:
        if len(specdir) == 0:
            specdir = os.path.join(SPECMATCHDIR, 'spectra')

        result.loc[:, 'snr'] = result.obs.apply(calc_snr, args=(specdir, ))
        result.sort_values(by='snr', ascending=False, inplace=True)

    return result


def calc_snr(obs, specdir):
    """Calculate the SNR of a given spectrum

    Args:
        obs (str): cps id of spectrum
        specdir (str): Directory to look for spectra
    Returns:
        float: 90th percentile SNR
    """
    try:
        f = fits.open(os.path.join(specdir, obs + '.fits'))
        serr = f[1].data
        snr = np.nanpercentile(1/serr, 90)
    except:
        print("Could not open " + obs)
        snr = np.nan

    return snr
