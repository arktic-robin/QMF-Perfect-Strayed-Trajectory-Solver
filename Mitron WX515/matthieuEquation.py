# This function was written with the purpose of being used for a Mechanical Engineering 3rd Year Dissertation at UoL(iverpool)
# The module code is ENGG301 titled "Individual Project"
# The author of the script is Bruno Reis
# Special thanks to my Project Supervisor Dave McIntosh for his expert advice

# Version: 515 -> Single Ion Calculations under Perfect Fields with Runge-Kutta 4
# Goals -> Accurately implement the Matthieu Equations and plot typical (XZ), (YZ) and (XY) spaces
# Supporting Function -> RK Solver
from math import cos


def matthieFunctionX(moveria, U, V, angularSpeed, time, x, startPhase=0):
    '''
    This is the X axis branch of the Matthieu Equations
    '''
    potentialWave = moveria * (0.5 * V * (1 + cos(angularSpeed * time + startPhase)) - U) * x
    return potentialWave

def matthieFunctionY(moveria, U, V, angularSpeed, time, y, startPhase=0):
    '''
    This is the X axis branch of the Matthieu Equations
    '''
    potentialWave = moveria * (U - 0.5 * V * (1 + cos(angularSpeed * time + startPhase))) * y
    print(cos(angularSpeed * time + startPhase))
    return potentialWave