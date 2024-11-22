# This is a program written for the ENGG309 Module of the Mechanical Engineering Programme at the University of Liverpool (Year 2024/25)
# The author of the program is called Bruno Reis
# Feel free to use it for whatever purpose you would like
import math
import random as rng
import numpy as np
import matplotlib.pyplot as plt
import qFieldPotential as qFP


# SET-UP
# Physical Parameter Definition
mass = 46
frequency = 1.3 * 1000
phaseSpeed = frequency * math.pi * 2
print(phaseSpeed)
potentialRatio = 0.95
radioPotential = 50
directPotential = radioPotential * potentialRatio
quadLength = 150
inscribedRadius = 6
xSrcOffset = 0
ySrcOffset = 0
axialSpeed = 5
initialPhase = 0
charge = 1
massRange = 1
spreadRange = 0.01
velocitySpread = 0.05

# Simulation Parameter Definition
phaseStepPactor = 1
phaseStep = phaseStepPactor * math.pi # Kinda redundant but at least I will be able to easily reference the precise value that I want
timeStep = 2 * phaseStep / phaseSpeed # Still very redundant
print(timeStep)
M = (frequency * quadLength) // (axialSpeed * phaseStepPactor)
M = int(M)
#M = (quadLength // axialSpeed) + 1

# Pre-allocate Memory Space
time = np.linspace(0, M * timeStep)
print(time.shape)
# use print(time) if you need to make sure that the pre-allocation is accurate

iPos = np.array([rng.uniform(-spreadRange, spreadRange) + xSrcOffset, rng.uniform(-spreadRange, spreadRange)  + ySrcOffset, 0])
iVel = np.array([rng.uniform(-velocitySpread, velocitySpread), rng.uniform(-velocitySpread, velocitySpread), axialSpeed]) # Another way by using the "constructor call"
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

#print(hstPos)
#print(hstVel)
#print(hstAcl)


# CALCULATIONS

# Trajectory Calculations

for i in range(M):
    #Loop Counting

    # Sum steps
    iAcl = qFP.qFP(mass, charge, inscribedRadius, radioPotential, directPotential, time[i], phaseSpeed) * iPos
    # use iAcl = np.array([ 1 , 2 , 0 ]) if you need to test the difference scheme
    iVel = iVel + iAcl * timeStep
    iPos = iPos + iVel * timeStep


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
