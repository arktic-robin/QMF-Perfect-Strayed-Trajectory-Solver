import math
import numpy as np
import qFieldPotential as qFP

def rk4(h, u, v, mass, charge, inscribedRadius, potDC, potRF, time, angularSpeed):
    # First cycle
    k1u = v
    k1v = qFP.qFP(mass, charge, inscribedRadius, potDC, potRF, time, angularSpeed)

    #Second Cycle
    k2u = v + k1v * h/2
    k2v = qFP.qFP(mass, charge, inscribedRadius, potDC, potRF, time + h/2, angularSpeed) * (u + k1u * h/2)

    # Third Cycle
    k3u = v + k2v * h/2
    k3v = qFP.qFP(mass, charge, inscribedRadius, potDC, potRF, time + h/2, angularSpeed) * (u + k2u * h/2)

    # Fourth Cycle
    k4u = v + k3v * h
    k4v = qFP.qFP(mass, charge, inscribedRadius, potDC, potRF, time + h/2, angularSpeed) * (u + k3u * h)

    # New Coordinates
    uNew = u + (k1u + 2 * (k2u + k3u) + k4u) * h/6
    vNew = v + (k1v + 2 * (k2v + k3v) + k4v) * h/6




    return uNew, vNew