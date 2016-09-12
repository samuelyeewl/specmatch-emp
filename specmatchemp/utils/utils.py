"""
@filename utils/utils.py

Various utility functions
"""

import numpy as np
from astropy import constants as c
from astropy import units as u


def calc_logg(radius, u_radius, mass, u_mass):
    """Calculates logg for a star from its mass and radius

    Args:
        radius: in Rsun
        u_radius: Uncertainty in radius
        mass: in Msun
        u_mass: Uncertainty in mass

    Returns:
        logg: in CGS
        u_logg: Propagated uncertainty
    """
    rstar = radius * c.R_sun
    mstar = mass * c.M_sun
    logg = np.log10((c.G * mstar/(rstar**2)).cgs.value)
    u_logg = (u_mass/mass + 2*u_radius/radius)*logg

    return logg, u_logg


def calc_radius(plx, u_plx, theta, u_theta):
    """Calculates stellar radius from parallax and angular diameter

    Args:
        plx: Parallax in mas
        u_plx: Uncertainty in parallax
        theta: Angular diameter in mas
        u_theta: Uncertainty in angular diameter

    Returns:
        radius: in Rsun
        u_radius: Uncertainty in radius
    """
    # convert units
    dist = (plx*u.marcsec).to(u.m, equivalencies=u.parallax())
    dimless_theta = (theta*u.marcsec).to('', equivalencies=u.dimensionless_angles())

    radius = (dist * dimless_theta / 2 / c.R_sun).value
    u_radius = (u_plx/plx + u_theta/theta) * radius

    return radius, u_radius

def calc_residuals(s1, w1, s2, w2):
    """Find the residuals between two spectra (s1-s2) when they are on the same
    wavelength scale but not necessarily within the same range

    Args:
        s1, w1: Spectrum 1
        s2, w2: Spectrum 2
    Returns:
        res, resw: Residuals
    """
    wmin = max(w1[0], w2[0])
    wmax = min(w1[-1], w2[-1])
    mask1 = np.logical_and(w1 > wmin, w1 < wmax)
    s1_masked = s1[mask1]
    mask2 = np.logical_and(w2 > wmin, w2 < wmax)
    s2_masked = s2[mask2]
    w_masked = w2[mask2]

    resid = s1_masked - s2_masked

    return resid, w_masked


def extend_array(arr, length, fill=np.nan):
    """Extends a numpy array to the given length, filling with provided fill
    value.

    Args:
        arr (np.ndarray): Array
        length (int): New length to extend to
        fill (float, optional): Fill value
    Returns:
        np.ndarray
    """
    new_arr = np.empty(length)
    if length < len(arr):
        new_arr[0:length] = arr[0:length]
    else:
        new_arr[0:len(arr)] = arr
        new_arr[len(arr):length] = fill

    return new_arr
