import numpy as np
import csv
import math

class QMF:
    """
    This class is intended to store the parameters of the filter
    """
    
    def __init__(self, name, tag, param=[]):
        """
        Initialise the filter object with a csv file
        or through the manual definition of the QMF
        parameters.

        Parameters
        ----------

        name: string
            The csv file to read the parameters
            from.
        tag: string
            The filter to be used from the QMF
            set in the csv file.
        param: list
            {NOT OPTIMAL} Manual definition of
            the parameters of the quadrupole.
            The name and tag must be set to
            None.

        Return
        ----------

        None
        """
        
        # Input Parameters
        self.name = name
        self.param = None
        self.tag = tag
        
        # Manual Parameter definition 
        if name == None:
            self.tag = param[0]
            self.radioPotential = param[1]
            self.directPotential = param[2]
            self.pureFrequency = param[3] * 1e6
            self.iniPhase = param[4]
            self.inscRadius = param[5] * 1e-3
            self.span = param[6] * 1e-3
            self.entryGap = param[7] * 1e-3
            self.exitGap = param[8] * 1e-3
            self.injSpeed = param[9] * 1e3

            self.angFrequency = 2 * math.pi * self.pureFrequency
            self.refVelocity = 2 * math.pi * self.pureFrequency * self.inscRadius
            self.hmc() # Determine the Hunter-McIntosh constants for this filter
    
    def hmc(self):
        """
        Calculate the Hunter-McIntosh fringing field
        constants based on a linear curve for the 1st
        coefficient and a piecewise function for
        the 2nd coefficient; this is done twice for
        the entry and exit fringing fields regions.

        Parameters
        ----------

        None

        Return
        ----------

        None
        """

        # Determine the enntry gap in units
        # of r_0 (dimensionless number)
        d1 = self.entryGap / self.inscRadius
        d2 = self.exitGap / self.inscRadius

        # 1st coefficients -> Linear Model
        self.a1 = 2.6688 - 2.3383 * d1
        self.a2 = 2.6688 - 2.3383 * d2

        
        # 2nd coefficient -> Piecewise

        if d1 < 0.125:
            self.b1 = 1.47
        elif d1 > 0.125 and d1 <= 0.25:
            self.b1 = 1.47 + (d1 - 0.125) * (8/125)
        elif d1 > 0.25 and  d1 <= 0.5:
            self.b1 = 1.55 + (d1 - 0.25) * (4/100)
        elif d1 > 0.5 and d1 <= 1:
            self.b1 = 1.56 - (d1 - 0.5) * (102/100)
        else:
            self.b1 = 1
        
        if d2 < 0.125:
            self.b2 = 1.47
        elif d2 > 0.125 and d2 <= 0.25:
            self.b2 = 1.47 + (d2 - 0.125) * (8/125)
        elif d2 > 0.25 and  d2 <= 0.5:
            self.b2 = 1.55 + (d2 - 0.25) * (4/100)
        elif d2 > 0.5 and d2 <= 1:
            self.b2 = 1.56 - (d2 - 0.5) * (102/100)
        else:
            self.b2 = 1
    
    def load_csv(self):
        """
        Read the filter set and reassign the values
        of the filter attributes.

        Parameters
        ----------

        None

        Return
        ----------

        None
        """
        
        with open("filter/" + self.name) as f:
            rows = csv.reader(f, delimiter=",")
            for row in rows:
                if row[0] == self.tag:
                    self.radioPotential = float(row[1])
                    self.directPotential = float(row[2])
                    self.pureFrequency = float(row[3]) * 1e6
                    self.iniPhase = float(row[4])
                    self.inscRadius = float(row[5]) * 1e-3
                    self.span = float(row[6]) * 1e-3
                    self.entryGap = float(row[7]) * 1e-3
                    self.exitGap = float(row[8]) * 1e-3
                    self.injSpeed = float(row[9])* 1e3

                    self.angFrequency = 2 * math.pi * self.pureFrequency
                    self.refVelocity = 2 * math.pi * self.pureFrequency * self.inscRadius
                    self.hmc()