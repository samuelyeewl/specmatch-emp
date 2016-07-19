"""
@filename analysis.py

Helper functions for analysis of results
"""

import numpy as np
import pandas as pd

from specmatchemp import library

def generate_sm_values(lib, results, method='lincomb', suffix='_sm'):
    """Generate the derived values and add it to the parameters table

    Args:
        lib (library.Library): Library object
        results (pd.DataFrame): results table 
        method ('lincomb', 'best_match', 'average'): Specify the method used
            - lincomb uses the weighted average of the parameters of a list of spectra
                where the coefficients were generated from the MatchLincomb procedure
                results table must contain ref_idxs, coeffs columns
            - best_match uses the parameters of the best matched spectrum
                results table must contain ref_idx column
            - average uses the unweighted average of the best matched spectra
                number of spectra to use can be specified by num_mean
        suffix (str): suffix to append to column name

    Returns:
        params (pd.DataFrame): parameter table
    """
    params = lib.library_params

    # Use linear combination results as sm values
    if method == 'lincomb':
        results.set_index('targ_idx', inplace=True)
        for p in library.STAR_PROPS:
            psm = p+suffix
            params[psm] = params.lib_index.apply(lambda i: \
                lincomb_props(lib.library_params, p, results.loc[i].ref_idxs, results.loc[i].coeffs))

    elif method == 'best_match':
        grouped_results = results.groupby('targ_idx')
        params['best_match'] = params.lib_index.apply(\
            lambda i: grouped_results.get_group(i).sort_values(by='chi_squared').iloc[0].ref_idx)
        for p in library.STAR_PROPS:
            psm = p+suffix
            params[psm] = params.best_match.apply(lambda i: params.loc[i, p])
        
        params['best_chi_squared'] = params.lib_index.apply(\
            lambda i: grouped_results.get_group(i).sort_values(by='chi_squared').iloc[0].chi_squared)

    elif method == 'average':
        grouped_results = results.groupby('targ_idx')
        params['best_n'] = params.lib_index.apply(lambda i: \
            np.array(grouped_results.get_group(i).sort_values(by='chi_squared').iloc[0:num_mean].ref_idx))
        for p in library.STAR_PROPS:
            psm = p+suffix
            params[psm] = params.best_n.apply(lambda i: params.loc[i, p].mean())
    
    return params


def lincomb_props(library_params, prop, idxs, coeffs):
    """Generates the weighted average of a given property

    Args:
        library_params (pd.DataFrame): Parameter table
        prop (str): Name of property column
        idxs (np.array): List of indices of reference spectra
        coeffs (np.array): Coefficients for weighted average
    Returns:
        sm_prop (np.float): Weighted average
    """
    assert np.isclose(np.sum(coeffs), 1, rtol = 1e-3, atol=1e-3), 'Coefficients must sum to 1'
    assert np.all(np.isfinite(coeffs)), 'Coefficients must be finite'
    assert np.all(np.isfinite(library_params.loc[idxs, prop])), 'Parameters must be finite'

    sm_prop = 0
    for i in range(len(idxs)):
        lib_prop = library_params.loc[idxs[i], prop]
        sm_prop += lib_prop*coeffs[i]
    return sm_prop

def generate_residuals(lib, suffix='_sm'):
    """Calculates the residuals between the derived and true values

    Args:
        lib (library.Library): Library object
        suffix (str): suffix to append to column name to obtain derived value

    Returns:
        params (pd.DataFrame): parameter table
    """
    params = lib.library_params

    for p in library.STAR_PROPS:
        presid = p+suffix+'_resid'
        psm = p+suffix
        lib.library_params[presid] = lib.library_params[psm] - lib.library_params[p]

    return params

def dist(star1, star2):
    """Distance between two stars in parameter space
    Normalized by Teff/100K, logg/0.1 dex, feh/0.1 dex
    
    Args:
        star1 (pd.DataFrame): Row of star 1
        star2 (pd.DataFrame): Row of star 2
    Returns:
        dist**2: Square of distance between the two stars
    """
    diff_teff = ((star1.Teff - star2.Teff)/100)**2
    diff_logg = ((star1.logg - star2.logg)/0.1)**2
    diff_feh = ((star1.feh - star2.feh)/0.1)**2
    return diff_teff + diff_logg + diff_feh


def find_closest_star(row, lib):
    """Helper function to find the closest star to the given star
    """
    return lib.library_params.apply(dist, args=(row,), axis=1).sort_values().index[1]

