import math
import numpy as np
from numpy import ones, zeros
import random as rng
import csv

import time as tm

class Gizmo:

    def __init__(self):

        # Keep this here for now, I might work on this later
        self.NUMBER_OF_MODELS = 3 # this doesn't brick the program
        self.BOOM_BOOST = 100000000000

        # ID purposes
        self.tag = None
        self.family = None

        # Storing crucial stuff here
        self._dst = []
        self._val = []

        # Parameters
        self.radioPotential = None
        self.directPotential = None
        self.pureFrequency = None
        self.angFrequency = None
        self.inscRadius = None
        self.maxPotential = None
    
    def hmc(self):

        # Determine the enntry gap in units
        # of r_0 (dimensionless number)
        d1 = 20 / 60
        d2 = 20 / 60

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
    
    def build(self, state, N):
        
        self.family = state
        
        # Useful, universal objects
        self.m = np.arange(1, self.NUMBER_OF_MODELS + 2)
        self.M = self.m[:, np.newaxis]
        self.M = np.repeat(self.M, N, axis=1)

        self.chargerTensor = np.empty((self.NUMBER_OF_MODELS + 1, 3, N))
        self.ZERO_MATRIX = np.zeros([3, N])
    
    def load(self, name):
    
        with open(name) as f:
            
            rows = csv.reader(f, delimiter=",")
            
            for row in rows:
                
                if row[0] == self.tag:
                    
                    self.radioPotential = float(row[1])
                    self.directPotential = float(row[2])
                    self.pureFrequency = float(row[3]) * 1e6
                    self.inscRadius = float(row[4]) * 1e-3
                    self.maxPotential = float(row[5])

                    self.angFrequency = 2 * math.pi * self.pureFrequency
                    self.hmc()

    def period(self):
        
        return 1 / self.pureFrequency
    
    def radialCheck(self, x, y):
        
        rad = self.inscRadius
        l = (abs(x) > rad) | (abs(y) > rad)
        
        return ~l
    
    def zoneCheck(self, z):

        D = self._dst
        
        # D is already in SI Units
        a = D[0]
        b = D[1]
        c = D[2]

        #print(a, b, c)

        l_1 = (z >= 0) & (z < a)

        l_2 = (z >= a) & (z < b)

        l_3 = (z >= b) & (z < c)

        l_4 = (z >= c)

        l = l_1 + (l_2 * 2) + (l_3 * 3) + (l_4 * 4)

        return l
    
    def pulse(self, x, time):

        state = self.family
        
        # Extract the l-vector (of shape 1xN)
        l = state.state[0, : ]


        # Determine the verification matrix (if zero, model stays)
        D = self.M - l
        V = 1 / (1 + self.BOOM_BOOST * D**2)

        # Determine the  force cards for all three (four) models {1, 2, 3, 4}
        F_1 = self._val[0].compute(x, time)
        F_2 = self._val[1].compute(x, time)
        F_3 = self._val[2].compute(x, time)
        F_4 = F_1 * 0
        
        # Assemble this pulse instance charger tensor
        self.chargerTensor[0, :, :] = F_1
        self.chargerTensor[1, :, :] = F_2
        self.chargerTensor[2, :, :] = F_3
        self.chargerTensor[3, :, :] = F_4

        # Determine the Transformation Matrix
        B = np.einsum('mij,mj->ij', self.chargerTensor, V)


        return B


class Generic_Zone:

    def __init__(self, master, span):

        self.master = master
        
        self.span = span * 1e-3
    
    def compute(self, state):
        return state * 0

class Detector_Zone(Generic_Zone):
    
    def __init__(self, master):
        self.master = master
        
        self.span = 0
    
    def compute(self, x, time):
        
        return self.master.ZERO_MATRIX

class Ideal_Zone(Generic_Zone):

    def __init__(self, master, span):
        
        super().__init__(master, span)

    def compute(self, x, time):

        # Reference to state object
        state = self.master.family

        X = x[0, :]
        Y = x[1, :]
        Z = x[2, :]

        #print(Z)

        ceta = state.iniPhase
        q = state.charge
        m = state.mass

        U = self.master.directPotential
        R = self.master.radioPotential
        w = self.master.angFrequency
        r = self.master.inscRadius

        e = q / m
        J = (U - R * math.cos(w * time + ceta)) / (r**2)

        A_X = (-1) * e * J * X
        A_Y = (1) * e * J * Y
        A_Z = (0) * e * J * Z

        y = np.vstack([A_X, A_Y, A_Z])
        
        return y

class Dawson_Entry_Zone(Generic_Zone):

    def __init__(self, master, span):
        
        super().__init__(master, span)

    def compute(self, x, time):

        # Zone definitions
        d = self.span

        # Reference to state object
        state = self.master.family

        X = x[0, :]
        Y = x[1, :]
        Z = x[2, :]

        ceta = state.iniPhase
        q = state.charge
        m = state.mass

        U = self.master.directPotential
        R = self.master.radioPotential
        w = self.master.angFrequency
        r = self.master.inscRadius

        e = q / m
        J = (U - R * math.cos(w * time + ceta)) / (r**2)


        # Zone-specific mechanics
        a = Z / d
        B = ( (X**2) - (Y**2) ) / d


        A_X = (-1) * e * J * np.multiply(X, a) 
        A_Y = (1) * e * J * np.multiply(Y, a)
        A_Z = (-1/2) * e * J * (B)

        y = np.vstack([A_X, A_Y, A_Z])

        return y



class Dawson_Exit_Zone(Generic_Zone):

    def __init__(self, master, span):
        
        super().__init__(master, span)

    def compute(self, x, time):

        # Zone definitions
        d = self.span

        # Reference to state object
        state = self.master.family
        l = self.master._dst[-2]

        X = x[0, :]
        Y = x[1, :]
        Z = x[2, :] - l # ensures that x3 starts at 0

        #print(Z, l)

        ceta = state.iniPhase
        q = state.charge
        m = state.mass

        U = self.master.directPotential
        R = self.master.radioPotential
        w = self.master.angFrequency
        r = self.master.inscRadius

        e = q / m
        J = (U - R * math.cos(w * time + ceta)) / (r**2)


        # Zone-specific mechanics
        a = 1 - (Z / d)
        B = ( (X**2) - (Y**2) ) / d


        A_X = (-1) * e * J * np.multiply(X, a) 
        A_Y = (1) * e * J * np.multiply(Y, a)
        A_Z = (1/2) * e * J * (B)

        y = np.vstack([A_X, A_Y, A_Z])

        return y





class HM_Entry_Zone(Generic_Zone):

    def __init__(self, master, span):
        
        super().__init__(master, span)

    def compute(self, x_vector, time):

        def f(x):

            z = x / self.span
            a = self.master.a1
            b = self.master.b1

            power = 0 - (a * z + b * (z**2))
            return 1 - np.exp(power)
    
        def f_prime(x):

            g = self.span
            z = x / g
            a = self.master.a1
            b = self.master.b1

            t1 = (a / g) + 2 * (b / (g**2)) * x
            t2 = np.exp(0 - (a * z + b * (z**2)))

            return np.multiply(t1, t2)

        # Zone definitions
        d = self.span

        # Reference to state object
        state = self.master.family

        # State definitions
        X = x_vector[0, :]
        Y = x_vector[1, :]
        Z = x_vector[2, :]

        ceta = state.iniPhase
        q = state.charge
        m = state.mass

        # Gizmo definitions
        U = self.master.directPotential
        R = self.master.radioPotential
        w = self.master.angFrequency
        r = self.master.inscRadius

        # Universal definitions
        e = q / m
        J = (U - R * math.cos(w * time + ceta)) / (r**2)


        # Zone-specific mechanics
        Alpha = f(Z)
        Beta = f_prime(Z)
        Delta = ( (X**2) - (Y**2) ) / d

        A_X = (-1) * e * J * np.multiply(X, Alpha)
        A_Y = (1) * e * J * np.multiply(Y, Alpha)
        A_Z = (-1/2) * e * J * np.multiply(Delta, Beta)


        # Output Prep
        y = np.vstack([A_X, A_Y, A_Z])
        
        return y

class HM_Exit_Zone(Generic_Zone):

    def __init__(self, master, span):
        
        super().__init__(master, span)

    def compute(self, x_vector, time):

        def f(x):

            z = x / self.span
            a = self.master.a1
            b = self.master.b1

            power = 0 - (a * z + b * (z**2))
            return 1 - np.exp(power)
    
        def f_prime(x):

            g = self.span
            z = x / g
            a = self.master.a1
            b = self.master.b1

            t1 = (a / g) + 2 * (b / (g**2)) * x
            t2 = np.exp(0 - (a * z + b * (z**2)))

            return np.multiply(t1, t2)

        # Zone definitions
        d = self.span

        # Reference to state object
        state = self.master.family

        # State definitions
        X = x_vector[0, :]
        Y = x_vector[1, :]
        Z = x_vector[2, :]

        ceta = state.iniPhase
        q = state.charge
        m = state.mass

        # Gizmo definitions
        U = self.master.directPotential
        R = self.master.radioPotential
        w = self.master.angFrequency
        r = self.master.inscRadius

        # Universal definitions
        e = q / m
        J = (U - R * math.cos(w * time + ceta)) / (r**2)


        # Zone-specific mechanics
        Alpha = f(Z)
        Beta = f_prime(Z)
        Delta = ( (X**2) - (Y**2) ) / d

        A_X = (-1) * e * J * np.multiply(X, Alpha)
        A_Y = (1) * e * J * np.multiply(Y, Alpha)
        A_Z = (1/2) * e * J * np.multiply(Delta, Beta)


        # Output Prep
        y = np.vstack([A_X, A_Y, A_Z])
        
        return y