import numpy as np
from numpy import fft
from scipy.special import jn
from scipy import integrate

def rotmacro_ft(sigma, xi, vsini, u1=0.43, u2=0.31, intres=100):
    """
    Fourier Transform of Rotational-Macroturbulent Broadening Kernel

    Hirano et al. (2011) equation B12, defines the Fourier Transform of the
    combined ratational and macroturbulent kernel.
    
    Args:
        sigma (float) : velocity frequency, s/km
        xi (float) : radial/tangential macroturblence (see Gray), km/s
        vsini (float) : projected surface velocity, km/s
        u1, u2 (float) : quadradtic limbdarkening parameters

    Returns
        Mt : Fourier transform of rotational-macroturbulent profile
    """
    
    tarr = np.outer(np.ones(sigma.shape[0]),np.linspace(0,1,intres)).transpose()


    t1 =( 
        (1 - u1 * ( 1 - np.sqrt( 1 - tarr**2 ) ) - 
         u2 * ( 1 - np.sqrt( 1 - tarr**2) )**2 ) / ( 1 - u1 / 3 - u2 / 6) 
    )

    t2 = (
        (
            np.exp(-np.pi**2 * xi**2 * sigma**2 * (1-tarr**2)) + 
            np.exp(-np.pi**2 * xi**2 * sigma**2 * tarr**2)
        ) * 
        jn(0,2*np.pi*sigma*vsini*tarr) * tarr
    )

    
    m = t1 * t2
    Mt = integrate.trapz(m,dx=np.mean(tarr[1:] - tarr[0:-1]),axis=0)
    return Mt

def rotmacro(n, dv, xi, vsini, **kwargs):
    """
    Rotational-Macroturbulent Broadening Kernel

    Computes the broadening profile for a range of velocity displacements.

    Args:
        n (int): number of points to evaluate kernel
        dv (float) : velocity spacing between points 
        xi (float) : macroturbulence
        vsini (float) : projected surface velocity, km/s
        **kwargs : passed to rotmacro_ft

    Returns
    -------
        varr (array): velocities shifts at which kernel is evaluated
        M (array): rotational macroturbulent kernel

    """
    nind = (n+1)/2 # number of independent points
    vmax = dv*nind
    varr = np.linspace(-vmax,vmax,n)

    # Fourier transform frequencies, for given velocity displacements
    sigarr = fft.fftfreq(n,d=dv)
    sigarr = sigarr[:nind]

    #M_ft_ep = np.array([rotmacro_ft(sig,xi,vsini,**kwargs) for sig in sigarr])
    M_ft = rotmacro_ft(sigarr,xi,vsini,**kwargs)

    M_ft = np.hstack([ M_ft, M_ft[1:][::-1] ])
    M = fft.ifft(M_ft)
    M = fft.fftshift(M)
    M = M/np.sum(M)
    return varr,M
