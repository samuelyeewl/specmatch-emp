"""
@filename match.py

Defines the Match class
"""
import pandas as pd



# lmfit python package
#

class Match:
    def calculate_chisquared(self, params):
        """
        Calculate the reduced chi-squared value
        """
        specmerged = pd.merge(self.target, self.reference, how='inner', on='w')
        residualsq = (specmerged['s_x']-specmerged['s_y'])**2
        ivar = 1/(specmerged['serr_x']**2+specmerged['serr_y']**2)
        return (residualsq*ivar).sum()/residualsq.count()

    def two_spectra(self, params):
        data = # observed spectrum
        tweaked_model = model * function(params)  = # library spectrum
        return data, tweaked_model

    def residual(self, params):
         data, tweaked_model =  two_spectra(params)
         # code happens
         return _residual




    def __init__(self, target, reference):
        """
        Constructor for the Match class

        target, reference spectra should be given as Pandas dataframes.
        """
        self.target = target
        self.reference = reference
        self.calculate_chisquared()
