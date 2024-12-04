# This program was written with the purpose of being used for a Mechanical Engineering 3rd Year Dissertation at UoL(iverpool)
# The module code is ENGG301 titled "Individual Project"
# The author of the script is Bruno Reis
# Special thanks to my Project Supervisor Dave McIntosh for his expert advice

# Version: 516 -> Single Ion Calculations under Perfect Fields with Runge-Kutta 4
# Goals -> Accurately implement the Matthieu Equations and plot typical (XZ), (YZ) and (XY) spaces

version = 517

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

# PHYSICAL PARAMETERS

# Welcome Message
welcomeMsg = f"""
This is a Quadrupole Mass Filter Trajectory Calculator, it calculates the numerical solution
to the Matthieu Equations of a Perfect Filter, This is version {version}. You will shortly be
prompted for a set of parameters for the simulation.
"""
print(welcomeMsg)

simParameters = {"Number of Ions" : ["Dimensionless", 1],
                 "Ion Mass" : ["Dalton", 46],
                 "Mass Spread" : ["Dalton", 0],
                 "Charge" : ["Multiple of Electron Charge", 1],
                 "Frequency" : ["Megahertz", 1.3],
                 "Radio Potential": ["Volt", 40452],
                 "Direct Potential" : ["Volt", 6552],
                 "Inscribed Radius" : ["Millimeter", 60],
                 "Quadrupole Span" : ["Millimeter", 150],
                 "Runge-Kutta Step Factor" : ["Dimensionless", 0.01],
                 "Ion Axial Speed" : ["kilometers-per-second", 4.52],
                 "Initial Phase" : ["radian", 0],
                 "Horizontal Displacement" : ["Millimeter", 2],
                 "Horizontal Spread" : ["Millimeter", 0],
                 "Vertical Displacement" : ["Millimeter", 2],
                 "Vertical Spread" : ["Millimeter", 0]}

# Loop through each parameter in the Dict
for p in simParameters:
  isChosen = False
  # Condition in case, the user enters a wrongly formatted string
  while (not isChosen):
    print(f"In terms of {p} in {simParameters[p][0]} units, what value should it be?")
    ans = input("Enter value: ")
    try:
       if ans == "":
          isChosen = True
          ans = simParameters[p][1]
       else:
          ans = float(ans)
          simParameters[p][1] = ans
          isChosen = True
    except:
       print("Uh-Oh, How unfortunate, Uh-Oh, How unfortunate!")
       print("Try entering the value again!")


# Parameter Assignment and Unit Conversion -> Check the Variable Audit
# Primary
N = simParameters["Number of Ions"][1]
mass = simParameters["Ion Mass"][1] * 1.66e-27
massSpread = simParameters["Mass Spread"][1] * 1.66e-27
charge = simParameters["Charge"][1] * 1.6e-19
freq = simParameters["Frequency"][1] * 1e6
potRF = simParameters["Radio Potential"][1]
potDC = simParameters["Direct Potential"][1]
insRad = simParameters["Inscribed Radius"][1] * 1e-3
span = simParameters["Quadrupole Span"][1] * 1e-3
hactor = simParameters["Runge-Kutta Step Factor"][1]
axSpeed = simParameters["Ion Axial Speed"][1] * 1e3
startPhase = simParameters["Initial Phase"][1]
xDispl = simParameters["Horizontal Displacement"][1] * 1e-3
xSpread = simParameters["Horizontal Spread"][1] * 1e-3
yDispl = simParameters["Vertical Displacement"][1] * 1e-3
ySpread = simParameters["Vertical Spread"][1] * 1e-3



# Secondary
angVel = 2 * math.pi * freq
h = hactor / freq
print(xDispl, yDispl, mass, h, charge, axSpeed, startPhase)

# Stability Numbers
a = 2 * potDC * charge / (mass * ( (insRad * angVel) ** 2))
q = 4 * potRF * charge / (mass * ( (insRad * angVel) ** 2))


# Array Declaration
iPos = np.empty( (3,), dtype=np.float64)
iVel = np.empty( (3,), dtype=np.float64)

hstPos = np.empty( (3, 3000), dtype=np.float64)
hstVel = np.empty( (3, 3000), dtype=np.float64)

hstTime = np.zeros( (3000,), dtype=np.float64)


# Initial Conditions -> Can link this to the terminal later
iPos[0] = xDispl # X initial Position
iPos[1] = yDispl # Y initial Position
iPos[2] = 0e-3 # Z initial Position

iVel[0] = 0e-3 # X initial Velocity
iVel[1] = 0e-3 # Y initial Velocity
iVel[2] = axSpeed # Z initial Velocity

# Setting the current coordinates to its appropriate index in the history array
hstPos[ :, 0] = iPos
hstVel[ :, 0] = iVel


"""
Computations
"""

period = 0
for i in range(2999):
    # Time & Pulse Wave Calculations
    period = i * h


    # Conditionals
    # Not implemented yet

    # Inverse Inertia
    moveria = charge / ( mass * ( insRad ** 2 ) )

    
    # Runge-Kutta 4  Solver of a Second Order Differential Equation + Euler Z-axis Solver
    x, y, u, v = rk4(i, moveria, iPos, iVel, h, period, potDC, potRF, angVel, startPhase)
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
ax1.plot(hstTime, x1)
ax1.set_title("XZ Space")
ax2.plot(hstTime, x2)
ax2.set_title("YZ Space")

plt.tight_layout()
plt.show()
