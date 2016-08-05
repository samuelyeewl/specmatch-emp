"""
@filename spectrum.py

Defines a spectrum object for use in SpecMatch-Emp
"""

import numpy as np 

class Spectrum():
    """Spectrum class

    This object is a container for a spectrum and other
    properties which may be required.
    """
    def __init__(self, w, s, serr=None, mask=None, name="", attrs={}):
        """
        Args:
            w (np.ndarray): Wavelength scale
            s (np.ndarray): Spectrum
            serr (np.ndarray): (optional) Error in spectrum
            mask (np.ndarray): (optional) Boolean array to mask out telluric lines
            name (str): (optional) Name associated with spectrum
            attrs (dict): (optional) Any further attributes
        """
        self.w = w
        self.s = s
        self.serr = serr
        self.mask = mask
        self.name = name
        self.attrs = attrs

    