"""Module for different broadening kernels.
"""
from scipy.special import jn
from scipy import integrate
import numpy as np
from numpy import fft

def rotmacro_ft(sigma, xi, vsini, u1=0.43, u2=0.31, intres=100):
    """Fourier Transform of Rotational-Macroturbulent Broadening Kernel

    Forumula B12 from Hirano et al. (2011)

    Args:
        sig (float): velocity frequency, s/km
        xi (float): radial/tangential macroturblence (see Gray), km/s
        vsini (float): projected surface velocity, km/s
        u1,u2 (float): quadradtic limbdarkening parameters

    Returns:
        array: Fourier transform of rot-macro profile at sig
    """ 
    
    tarr = np.outer( np.ones(sigma.shape[0]), np.linspace(0, 1, intres) )
    tarr = tarr.transpose()

    t1 = (
        (1 - u1 * ( 1 - np.sqrt(1 - tarr**2) ) - 
         u2 * ( 1 - np.sqrt( 1 - tarr**2 ) )**2) / ( 1 - u1 / 3 - u2 / 6 )
    )
    t2 = (
        ( np.exp( -np.pi**2 * xi**2 * sigma**2 * (1-tarr**2) ) + 
          np.exp( -np.pi**2 * xi**2 * sigma**2 * tarr**2) ) * 
        jn(0,2*np.pi*sigma*vsini*tarr) * tarr
    )

    m = t1 * t2
    kernel_ft = integrate.trapz(m,dx=np.mean(tarr[1:] - tarr[0:-1]),axis=0)
    return kernel_ft

def rotmacro(n, dv, xi, vsini, **kwargs):
    """Rotational-Macroturbulent Broadening Kernel

    Args:
        n (int): number of points over which to calculate kernel
        dv (float): spacing between points (km/s)
        xi (float): macroturbulence (km/s)
        vsini (float): projected rotational velocity
        **kwargs: passed to rotmacro_ft

    Returns:
        array: velocities
        array: broadening kernel (same size as varr)
    """
    nind = (n + 1) / 2 
    vmax = dv * nind
    varr = np.linspace( -vmax, vmax, n)

    # Fourier transform frequencies, for given velocity displacements
    sigarr = fft.fftfreq(n,d=dv)
    sigarr = sigarr[:nind]
    kernel_ft = rotmacro_ft(sigarr,xi,vsini,**kwargs)
    kernel_ft = np.hstack([ kernel_ft, kernel_ft[1:][::-1] ])
    kernel = fft.ifft(kernel_ft)
    kernel = fft.fftshift(kernel)
    kernel = kernel/np.sum(kernel)
    return varr,kernel

def rot(n, dv, vsini, u1=0.6):
    """Rotational Kernel

    Create a 1-d convolution kernel to broaden a spectrum from a rotating star
    
    Args:
        deltaV (float): spacing between kernel points (km/s)
        vsini (float): projecte rotational velocoity

        u1 (Optional[float]): linear limb-darkening
            coefficient. Default value of 0.6 is typical for
            photospheric lines. The specific intensity I at any angle
            theta from the specific intensity Icen at the center of
            the disk is given by:
  
            I = Icen*(1-epsilon*(1-cos(theta))

    Returns:
        array: The convolution kernel vector for the specified
            rotational velocity. kernel will always contain an odd number of 
            points

    Notes:
        Adapted from lsf_rotate.pro by Erik Petigura

    """
    nind = (n + 1) / 2 # number of independent points
    vmax = dv * nind
    varr = np.linspace(-vmax, vmax, n).astype(float)

    e1 = 2.0 * (1.0 - u1)
    e2 = np.pi * u1 / 2.0
    e3 = np.pi * (1.0 - u1 / 3.0)

    x = varr / vsini
    x1 = np.abs(1.0 - x**2)
    kernel = ( ( e1 * np.sqrt(x1) + e2 * x1 ) / e3).astype(float)
    kernel[np.abs(x) > 1] = 0 # kernel should be 0 for v > vsini 
    kernel /= kernel.sum() # normalize 
    return varr, kernel
