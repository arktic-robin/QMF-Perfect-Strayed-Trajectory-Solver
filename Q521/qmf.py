import numpy as np
import csv
import math

class QMF:
    """
    This class is intended to store the parameters of the filter
    """
    
    def __init__(self, name, tag, param=[]):
        # Input Parameters
        self.name = name
        self.param = None
        self.tag = tag
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
    
    def flowPeriod(self):
        # Important Calculations
        T1 = self.entryGap /  self.injSpeed
        T2 = self.span /  self.injSpeed
        T3 = self.exitGap /  self.injSpeed
        return [T1, T2, T3]
    
    def load_csv(self):
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