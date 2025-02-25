import numpy as np
import math

class Generic_Law:
    """
    This is a placeholder
    """

    def __init__(self, specialClass=None):
        self.special = specialClass
    
    def compute(self):
        # Also a placeholder, override this
        y = 131
        return y
    
class Matthieu_Equation(Generic_Law):
    
    def __init__(self, specialClass=None):
        super().__init__(specialClass)
    
    def compute(self, qmf, ion):
        x = ion.pos
        v = ion.vel

        c0 = np.zeros([3,])
        c0[0] = -1
        c0[1] = 1

        c1 = ion.charge / (ion.mass * (qmf.inscRadius ** 2))
        c2 = qmf.directPotential + qmf.radioPotential * math.sin(ion.phaseTime[1])


        x = np.transpose(x)
        y = c0 * c1 * c2 * x
        y = np.transpose(y)
        
        return y # Placeholder -> I will work on this later today or tomorrow