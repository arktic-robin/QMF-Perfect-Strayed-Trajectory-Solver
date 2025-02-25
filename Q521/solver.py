import numpy as np
import math

class Generic_Solver:
    """
    This is a template class for all the integration
    methods I will use, that is: RK4, Euler and RK4.5
    """

    def __init__(self, specialClass=None):
        self.special = specialClass

    def predict(self):
        # This is literally just for testing
        x = self.ion.pos + 1
        v = self.ion.vel + 2
        a = self.ion.acl + 9
        q = self.ion.gen + 12
        return x, v, a, q
    
class RK4(Generic_Solver):
    """
    This is the Runge-Kutta 4 Method
    """

    def __init__(self, specialClass=None):
        super().__init__(filter, specialClass)
    
    def predict(self):
        x = self.ion.pos + 1
        v = self.ion.vel + 2
        a = self.ion.acl + 9
        q = self.ion.gen + 12
        return x, v, a, q
    
class Euler(Generic_Solver):
    """
    This is the Runge-Kutta 4 Method
    """

    def __init__(self, specialClass=None):
        super().__init__(specialClass)
    
    def predict(self, law, h, qmf, ion):
        a = law.compute(qmf, ion)
        v = ion.vel + h * a
        x = ion.pos + h * v
        q = ion.gen + 12

        return x, v, a, q