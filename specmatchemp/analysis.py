"""
@filename analysis.py

Helper functions for analysis of results
"""

import numpy as np
import pandas as pd

from specmatchemp.library import Library

def generate_sm_values(params, results, method='lincomb', suffix='_sm', cscol='chi_squared', refcol='ref_idxs', coeffcol='coeffs'):
    """Generate the derived values and add it to the parameters table

    Args:
        params (pd.DataFrame): Paramter table
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

    # Use linear combination results as sm values
    if method == 'lincomb':
        results = results.set_index('targ_idx')
        for p in Library.STAR_PROPS:
            psm = p+suffix
            params.loc[:,psm] = params.lib_index.apply(lambda i: \
                lincomb_props(params, p, results.loc[i, refcol], results.loc[i, coeffcol]))

    elif method == 'best_match':
        grouped_results = results.groupby('targ_idx')
        params.loc[:,'best_match'+suffix] = params.lib_index.apply(\
            lambda i: grouped_results.get_group(i).sort_values(by=cscol).iloc[0].ref_idx)
        for p in Library.STAR_PROPS:
            psm = p+suffix
            params.loc[:,psm] = params['best_match'+suffix].apply(lambda i: params.loc[i, p])
        
        params.loc[:,'best_chi_squared'+suffix] = params.lib_index.apply(\
            lambda i: grouped_results.get_group(i).sort_values(by=cscol).iloc[0][cscol])

    elif method == 'average':
        grouped_results = results.groupby('targ_idx')
        params.loc[:,'best_n'+suffix] = params.lib_index.apply(lambda i: \
            np.array(grouped_results.get_group(i).sort_values(by=cscol).iloc[0:num_mean].ref_idx))
        for p in Library.STAR_PROPS:
            psm = p+suffix
            params.loc[:,psm] = params['best_n'+suffix].apply(lambda i: params.loc[i, p].mean())
    
    return params


def lincomb_props(params, prop, idxs, coeffs):
    """Generates the weighted average of a given property

    Args:
        params (pd.DataFrame): Parameter table
        prop (str): Name of property column
        idxs (np.array): List of indices of reference spectra
        coeffs (np.array): Coefficients for weighted average
    Returns:
        sm_prop (np.float): Weighted average
    """
    assert np.isclose(np.sum(coeffs), 1, rtol = 1e-3, atol=1e-3), 'Coefficients must sum to 1'
    assert np.all(np.isfinite(coeffs)), 'Coefficients must be finite'
    # assert np.all(np.isfinite(library_params.loc[idxs, prop])), 'Parameters must be finite'

    sm_prop = 0
    for i in range(len(idxs)):
        lib_prop = params.loc[idxs[i], prop]
        sm_prop += lib_prop*coeffs[i]
    return sm_prop

def generate_residuals(params, suffix='_sm', props=Library.STAR_PROPS):
    """Calculates the residuals between the derived and true values

    Args:
        params (pd.Dataframe): parameter table
        suffix (str): suffix of derived values

    Returns:
        params (pd.DataFrame): parameter table
    """
    for p in props:
        presid = p+suffix+'_resid'
        psm = p+suffix
        params.loc[:,presid] = params[psm] - params[p]

    # calculate delta r/r
    presid = 'dr_r'+suffix+'_resid'
    params.loc[:,presid] = (params['radius'+suffix]-params['radius'])/params['radius']

    return params

def detrend_params(params, suffix='_sm'):
    """Detrend the parameters

    Args:
        params (pd.Dataframe): parameter table
        suffix (str): suffix of derived values

    Returns:
        params (pd.DataFrame): parameter table
        polycoeffs (dict of ndarrays): Fit coeffs
            - 'Teff_cool': Teff fit for cool stars (Teff < 4500)
            - 'Teff_hot': Teff fit for hot stars (Teff > 4500)
            - 'feh': feh fit
    """
    polycoeffs = {}

    # Fit a linear trend to cool stars
    T_derived = 'Teff'+suffix
    T_detrend = 'Teff'+suffix+'_detrend'

    cool = params.query('Teff < 4500')
    p = np.polyfit(cool['Teff'], cool['Teff'+suffix+'_resid'], 1)
    # correction in coefficients for derived values
    p0 = p[0]/(p[0]+1)
    p1 = p[1]/(p[0]+1)
    params[T_detrend] = params[T_derived]
    params.loc[:,T_detrend] = params.apply(lambda row: \
        row[T_derived] - p0*row[T_derived] - p1 \
        if row[T_derived] < 4500*(p[0]+1) + p[1] else row[T_detrend], axis=1)
    polycoeffs['Teff_cool'] = p

    # Fit a separate linear trend to hot stars
    hot = params.query('Teff >= 4500')
    p = np.polyfit(hot['Teff'], hot['Teff'+suffix+'_resid'], 1)
    p0 = p[0]/(p[0]+1)
    p1 = p[1]/(p[0]+1)
    params.loc[:,T_detrend] = params.apply(lambda row: \
        row[T_derived] - p0*row[T_derived] - p1 \
        if row[T_derived] >= 4500*(p[0]+1) + p[1] else row[T_detrend], axis=1)
    polycoeffs['Teff_hot'] = p

    # Fit a trend to giant stars (R > 1.2 Rsun)
    giants = params.query('1.1 < radius < 2.5')
    giants = giants[np.logical_not(np.isnan(giants['radius'+suffix]))]
    p = np.polyfit(np.log(giants['radius']), giants['radius'+suffix+'_resid']/giants['radius'], 1)
    params.loc[:,'radius'+suffix+'_detrend'] = params.apply(lambda row:\
        row['radius'+suffix] - row['radius'+suffix]*(p[0]*np.log(row['radius'+suffix]) + p[1]) \
        if row['radius'+suffix] > 1.1 and row['radius'+suffix] < 2.5 else row['radius'+suffix], axis=1)
    polycoeffs['radius_giants'] = p

    # Fit a linear trend to feh
    p = np.polyfit(params['feh'], params['feh'+suffix+'_resid'], 1)
    p0 = p[0]/(p[0]+1)
    p1 = p[1]/(p[0]+1)
    params.loc[:,'feh'+suffix+'_detrend'] = params['feh'+suffix] - p0*params['feh'+suffix] - p1
    polycoeffs['feh'] = p

    params = generate_residuals(params, suffix+'_detrend', props=['Teff', 'radius', 'feh'])

    return (params, polycoeffs)


def dist(star1, star2):
    """Distance between two stars in parameter space
    Normalized by Teff/100K, DeltaR/R/0.2 dex, feh/0.1 dex
    
    Args:
        star1 (pd.DataFrame): Row of star 1
        star2 (pd.DataFrame): Row of star 2
    Returns:
        dist**2: Square of distance between the two stars
    """
    diff_teff = ((star1.Teff - star2.Teff)/100)**2
    # diff_logg = ((star1.logg - star2.logg)/0.1)**2
    diff_radius = ((star1.radius - star2.radius)/(np.average([star1.radius,star2.radius]))/0.2)**2
    diff_feh = ((star1.feh - star2.feh)/0.1)**2
    return diff_teff + diff_radius + diff_feh


def find_closest_star(row, lib):
    """Helper function to find the closest star to the given star
    """
    return lib.library_params.apply(dist, args=(row,), axis=1).sort_values().index[1]

