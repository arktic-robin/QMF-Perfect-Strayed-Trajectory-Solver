import numpy as np
import math
from qmf import QMF
from history import History
from law import Matthieu_Equation
from solver import Euler

class Corral:
    """This class is a composition and polymorphism heavy
    type which serves as a generic trajectory calculator
    for a given regime of the quadrupole system."""

    def __init__(self, solver, law, timeStep, steps, filter, history):
        self.h = timeStep
        self.steps = steps 
        self.qmf = filter
        self.ion = history

        self.solver = solver()
        self.law = law()
    
    def lockOn(self):
        for i in range(self.steps):
            i += 1
            
            self.ion.step(i, self.h, self.qmf.angFrequency)
            tempIon = self.ion.eject(i - 1)
            x, v, a, q = self.solver.predict(self.law, self.h, self.qmf, tempIon)

            self.ion.pos[:, i, :] += x
            self.ion.vel[:, i, :] += v
            self.ion.acl[:, i, :] += a
            self.ion.gen[:, i, :] += q

            #print(self.ion.pos[:, i, :])

            

    def main(self):
        t = self.ion.phaseTime
        x = self.ion.pos
        v = self.ion.vel
        a = self.ion.acl
        q = self.ion.gen
        return t, x, v, a, q