# This function was written with the purpose of being used for a Mechanical Engineering 3rd Year Dissertation at UoL(iverpool)
# The module code is ENGG301 titled "Individual Project"
# The author of the script is Bruno Reis
# Special thanks to my Project Supervisor Dave McIntosh for his expert advice

# Version: 515 -> Single Ion Calculations under Perfect Fields with Runge-Kutta 4
# Goals -> Accurately implement the Matthieu Equations and plot typical (XZ), (YZ) and (XY) spaces
# Supporting Function -> RK Solver

from matthieuEquation import matthieFunctionX, matthieFunctionY

def rk4(i, moveria, position, velocity, h, period, U, V, w, startPhase = 0):
    '''
    This function is an RK4 solver for the Mathieu Equation.

    Parameters:
    - moveria: Parameter for Mathieu function computations.
    - position: List or array containing initial positions [x, y, z].
    - velocity: List or array containing initial velocities [vx, vy, vz].
    - h: Time step size.
    - period: Current time period for the Mathieu equation.
    - U, V: Parameters for Mathieu function computations.
    - w: Angular frequency.
    - startPhase: Initial phase offset (default: 0).

    Returns:
    - u1New: Updated position in the X-axis.
    - u2New: Updated position in the Y-axis.
    - v1New: Updated velocity in the X-axis.
    - v2New: Updated velocity in the Y-axis.
    '''

    # Defining temporary variables
    u1 = position[0]
    u2 = position[1]
    u3 = position[2]

    v1 = velocity[0]
    v2 = velocity[1]
    v3 = velocity[2]





    # First Cycle
    # X - Axis
    k1u1 = v1
    k1v1 = matthieFunctionX(moveria, U, V, w, period, u1, startPhase)

    # Y - Axis
    k1u2 = v2
    k1v2 = matthieFunctionY(moveria, U, V, w, period, u2, startPhase)

    # Second Cycle
    # X - Axis
    k2u1 = v1 + k1v1 * h / 2
    k2v1 = matthieFunctionX(moveria, U, V, w, period + h / 2, u1 + k1u1 * h / 2, startPhase)

    # Y - Axis
    k2u2 = v2 + k1v2 * h / 2
    k2v2 = matthieFunctionY(moveria, U, V, w, period + h / 2, u2 + k1u2 * h / 2, startPhase)

    # Third Cycle
    # X - Axis
    k3u1 = v1 + k2v1 * h / 2
    k3v1 = matthieFunctionX(moveria, U, V, w, period + h / 2, u1 + k2u1 * h / 2, startPhase)

    # Y - Axis
    k3u2 = v2 + k2v2 * h / 2
    k3v2 = matthieFunctionY(moveria, U, V, w, period + h / 2, u2 + k2u2 * h / 2, startPhase)

    # Fourth Cycle
    # X - Axis
    k4u1 = v1 + k3v1 * h
    k4v1 = matthieFunctionX(moveria, U, V, w, period + h, u1 + k3u1 * h, startPhase)

    # Y - Axis
    k4u2 = v2 + k3v2 * h
    k4v2 = matthieFunctionY(moveria, U, V, w, period + h, u2 + k3u2 * h, startPhase)

    '''
    if i in [0, 1, 2]:  # Log every 100 iterations
        print(f"Iteration {i} | u, k1: {k1u2}, k2: {k2u2}, k3: {k3u2}, k4: {k4u2}")
        print(f"Iteration {i} | v, k1: {k1v2}, k2: {k2v2}, k3: {k3v2}, k4: {k4v2}")
        pass
    '''


    # Final Integration
    # X - Axis
    u1New = u1 + (k1u1 + (2 * (k2u1 + k3u1)) + k4u1) * h / 6
    v1New = v1 + (k1v1 + (2 * (k2v1 + k3v1)) + k4v1) * h / 6

    # Y - Axis
    u2New = u2 + (k1u2 + (2 * (k2u2 + k3u2)) + k4u2) * h / 6
    v2New = v2 + (k1v2 + (2 * (k2v2 + k3v2)) + k4v2) * h / 6



    return u1New, u2New, v1New, v2New