# This is a program written for the ENGG309 Module of the Mechanical Engineering Programme at the University of Liverpool (Year 2024/25)
# The author of the program is called Bruno Reis
# Feel free to use it for whatever purpose you would like
import math
import random as rng
import numpy as np
import matplotlib.pyplot as plt
import qFieldPotential as qFP
from rk4 import rk4


# SET-UP

#Constants
dalton2Kg = 1.6605 * (10**-27)

# Physical Parameter Definition
mass = 46
frequency = 1.3
phaseSpeed = frequency * math.pi * 2
print(phaseSpeed)
potentialRatio = 0.95
radioPotential = 25
directPotential = 12.5 # radioPotential * potentialRatio
quadLength = 150 / 1000
inscribedRadius = 6 / 1000
xSrcOffset = 0
ySrcOffset = 0
axialSpeed = 5 / 1000
initialPhase = 0
charge = 1
massRange = 1
spreadRange = 0.01
velocitySpread = 0.05

# Simulation Parameter Definition
timeStep = 0.005 # Still very redundant
print(timeStep)
M = (quadLength / axialSpeed) / timeStep
M = int(M)
period = 0
#M = (quadLength // axialSpeed) + 1

# Pre-allocate Memory Space
# use print(time) if you need to make sure that the pre-allocation is accurate

#iPos = np.array([rng.uniform(-spreadRange, spreadRange) + xSrcOffset, rng.uniform(-spreadRange, spreadRange)  + ySrcOffset, 0])
#iVel = np.array([rng.uniform(-velocitySpread, velocitySpread), rng.uniform(-velocitySpread, velocitySpread), axialSpeed]) # Another way by using the "constructor call"

iPos = np.array([0.0045, 0.0045, 0])
iVel = np.array([0.001, 0.001, axialSpeed])
iAcl = np.zeros([1, 3]) # One way of making arrays

hstPos = np.zeros([M, 3])
hstVel = np.zeros([M, 3])
hstAcl = np.zeros([M, 3])

hstPos[0, :] = iPos
hstVel[0, :] = iVel
hstAcl[0, :] = iAcl

print(M)
print(iPos)
print(iVel)
print(iAcl)


# CALCULATIONS

# Trajectory Calculations

for i in range(M):
    #Loop Counting
    period += timeStep

    # Solver Steps
    iPos, iVel = rk4(timeStep, iPos, iVel, mass, charge, inscribedRadius, radioPotential, directPotential, period, phaseSpeed)



    # Annotate Coordinates
    hstPos[ i , : ] = iPos
    hstVel[ i , : ] = iVel
    hstAcl[ i , : ] = iAcl



# DATA VISUALIZATION
# Pull Path Columns

x1 = hstPos[ : , 0 ]
x2 = hstPos[ : , 1 ]
x3 = hstPos[ : , 2 ]

v1 = hstVel[ : , 0 ]
v2 = hstVel[ : , 1 ]
v3 = hstVel[ : , 2 ]

a1 = hstAcl[ : , 0 ]
a2 = hstAcl[ : , 1 ]
a3 = hstAcl[ : , 2 ]

# Create Figures

fig, (ax1, ax2) = plt.subplots(2, 1)
ax1.plot(x3, x1)
ax1.set_title("XZ Space")
ax2.plot(x3, x2)
ax2.set_title("YZ Space")

plt.tight_layout()
plt.show()
