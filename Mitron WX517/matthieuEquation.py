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
    This is the X-axis branch of the Mathieu Equations.

    Parameters:
    - moveria: Parameter scaling the potential wave.
    - U: Constant potential offset.
    - V: Amplitude of the periodic potential.
    - angularSpeed: Angular frequency of the periodic potential.
    - time: Time at which the potential is evaluated.
    - x: Position along the X-axis.
    - startPhase: Initial phase offset for the cosine wave (default: 0).

    Returns:
    - The potential wave contribution for the X-axis.
    '''
    timeDependentTerm = V * cos(angularSpeed * time + startPhase)
    potentialWave = moveria * (timeDependentTerm - U) * x
    #print(f" The Inertia: {moveria} | The Phase: {angularSpeed * time + startPhase} | The Pulse : {cos(angularSpeed * time + startPhase)} | The Wave: {timeDependentTerm - U} | The End: {x}")
    return potentialWave

def matthieFunctionY(moveria, U, V, angularSpeed, time, y, startPhase=0):
    '''
    This is the Y-axis branch of the Mathieu Equations.

    Parameters:
    - moveria: Parameter scaling the potential wave.
    - U: Constant potential offset.
    - V: Amplitude of the periodic potential.
    - angularSpeed: Angular frequency of the periodic potential.
    - time: Time at which the potential is evaluated.
    - y: Position along the Y-axis.
    - startPhase: Initial phase offset for the cosine wave (default: 0).

    Returns:
    - The potential wave contribution for the Y-axis.
    '''
    timeDependentTerm = V * cos(angularSpeed * time + startPhase)
    potentialWave = moveria * (U - timeDependentTerm) * y
    return potentialWave