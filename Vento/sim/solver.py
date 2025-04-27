import numpy as np
from numpy import ones, zeros
import random as rng
import csv

class Generic_Solver:

    def __init__(self, master):

        self.master = master

    def predict(self, gizmo, t):
        return 0

class Euler(Generic_Solver):

    def __init__(self, master):
        super().__init__(master)

    def predict(self, gizmo, t):

        master = self.master
        h = master.h

        # state variables
        X = master.state[1:4, :]
        V = master.state[4:7, :]

        A = gizmo.pulse(X, t)

        # Reassignment
        master.state[1:4, : ] += h * V
        master.state[4:7, : ] += h * A

class RK4(Generic_Solver):

    def __init__(self, master):
        super().__init__(master)
    
    def predict(self, gizmo, t):

        def acl(X, tim):

            # Runs pulse and gets the prediction 
            return gizmo.pulse(X, tim)

        master = self.master
        h = master.h

        # Initial estimate values 
        X = master.state[1:4, :]
        V = master.state[4:7, :]

        # RK4 for position and velocity
        k1v = acl(X, t)
        k1x = V

        k2v = acl(X + 0.5 * h * k1x, t + 0.5 * h)
        k2x = V + 0.5 * h * k1v

        k3v = acl(X + 0.5 * h * k2x, t + 0.5 * h)
        k3x = V + 0.5 * h * k2v

        k4v = acl(X + h * k3x, t + h)
        k4x = V + h * k3v


        # Reassignment
        master.state[1:4, : ] += (h/6) * (k1x + 2 * k2x + 2 * k3x + k4x)
        master.state[4:7, : ] += (h/6) * (k1v + 2 * k2v + 2 * k3v + k4v)