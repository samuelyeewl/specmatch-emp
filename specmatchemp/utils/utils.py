"""
@filename buildlib/utils.py

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
    u_logg = (row['e_M']/row['M'] + 2*row['e_R']/row['R'])*logg

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
    
