# This program was written with the purpose of being used for a Mechanical Engineering 3rd Year Dissertation at UoL(iverpool)
# The module code is ENGG301 titled "Individual Project"
# The author of the script is Bruno Reis
# Special thanks to my Project Supervisor Dave McIntosh for his expert advice

# Version: 515 -> Single Ion Calculations under Perfect Fields with Runge-Kutta 4
# Goals -> Accurately implement the Matthieu Equations and plot typical (XZ), (YZ) and (XY) spaces



"""
For the majority of the program, I will be using double-quoted docstrings to separate the file into logical chunks
This is not standard practice, but it makes it easier for me to work around all the logic and statements
Note that for Functions, I use single-quoted docstrings for their documentation
:)
"""

# Third-party Modules and Libraries
import math
import numpy as np
from matplotlib import pyplot as plt

# Coupled functions
from rk4 import rk4



"""
Pre-allocation and Definitions
"""
# Array Declaration
iPos = np.empty( (3,), dtype=np.float64) # Didn't know Numpy had float 128-bit, cool
iVel = np.empty( (3,), dtype=np.float64) # They even got 256-bit!!!

hstPos = np.empty( (3, 30000), dtype=np.float64)
hstVel = np.empty( (3, 30000), dtype=np.float64)

hstTime = np.zeros( (30000,), dtype=np.float64)

# Physical Parameters
h = 0.0005 # h is literally just delta t in seconds, in fact, everything is in SI units for 515
period = 0
insRad = 6e-3 # this is the inscribed radius of the filter, in meters
span = 1500e-3 # This is the length of the quadrupoles, also in meters
potDC = 20 # Volts
potRF = 40 # Volts again
startPhase = 0 # Phase in radians
freq = 1.3e15 # Hz, in future versions I will prompt the user for the Frequency in megahertz
angVel = 2 * math.pi * freq
axSpeed = 0.15e-3
mass = 2.301e-25 # The mass is in kg, the ion in question is Barium-138
charge = 2 * 1.9e-19 # Barium +2 ion

# Initial Conditions -> Can link this to the terminal later
iPos[0] = 2e-3 # X initial Position
iPos[1] = 2e-3 # Y initial Position
iPos[2] = 2e-3 # Z initial Position

iVel[0] = 0e-3 # X initial Velocity
iVel[1] = 0e-3 # Y initial Velocity
iVel[2] = axSpeed # Z initial Velocity

# Setting the current coordinates to its appropriate index in the history array
hstPos[ :, 0] = iPos
hstVel[ :, 0] = iVel


"""
Computations
"""

for i in range(round(span / iVel[2]) ):
    # Time & Pulse Wave Calculations
    #print(i * h)
    period = i * h
    pulse = math.cos(angVel * period + startPhase)


    # Conditionals
    """
    This can be blank for now, but when I implement a strayef/fringing fields routine, this is probably where it comes in
    """

    # Inverse Inertia
    #moveria = charge / ( mass * ( insRad ** 2 ) )
    moveria = 12

    
    # Runge-Kutta 4  Solver of a Second Order Differential Equation + Euler Z-axis Solver
    x, y, u, v = rk4(moveria, iPos, iVel, h, period, potDC, potRF, angVel, startPhase)
    z = period * iVel[2]

    iPos[0] = x
    iPos[1] = y
    iPos[2] = z

    iVel[0] = u
    iVel[1] = v
    

    hstPos[ :, i + 1] = iPos
    hstVel[ :, i + 1] = iVel
    hstTime[i + 1] = period

"""
Data Visualization
"""

x1 = hstPos[ 0 , : ]
x2 = hstPos[ 1 , : ]
x3 = hstPos[ 2 , : ]

v1 = hstVel[ 0 , : ]
v2 = hstVel[ 1 , : ]
v3 = hstVel[ 2 , : ]

# Create Figures

fig, (ax1, ax2) = plt.subplots(2, 1)
ax1.scatter(hstTime, x1)
ax1.set_title("XZ Space")
ax2.scatter(hstTime, x2)
ax2.set_title("YZ Space")

plt.tight_layout()
plt.show()
