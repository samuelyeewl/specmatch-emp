"""
@filename match.py

Defines the Match class
"""
import pandas as pd

class Match:
    def calculate_chisquared(self):
        """
        Calculate the reduced chi-squared value
        """
        specmerged = pd.merge(self.target, self.reference, how='inner', on='w')
        residualsq = (specmerged['s_x']-specmerged['s_y'])**2
        ivar = 1/(specmerged['serr_x']**2+specmerged['serr_y']**2)
        return (residualsq*ivar).sum()/residualsq.count()

    def __init__(self, target, reference):
        """
        Constructor for the Match class

        target, reference spectra should be given as Pandas dataframes.
        """
        self.target = target
        self.reference = reference
        self.calculate_chisquared()
