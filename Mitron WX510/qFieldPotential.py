import math
import numpy as np

def qFP(mass, charge, inscribedRadius, potDC, potRF, time, angularSpeed):
    t1 = charge / (mass * inscribedRadius ** 2) 
    t2 = potDC - potRF * math.cos(time * angularSpeed)
    np.array([ - t1*t2 , t1*t2 , 0 ])
    return np.array([ - t1*t2 , t1*t2 , 0 ])
