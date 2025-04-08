import numpy as np
import math

class Generic_Solver:
    """
    This is a template class for all the integration
    methods I will use, that is: RK4, Euler and RK2
    """

    def __init__(self, specialClass=None):
        self.special = specialClass

    def predict(self):
        """
        Placeholder for the solver methods.
        The 'predict' method should all have
        the exact same inputs.

        Parameters
        ----------

        law: Generic_Law (& sub-classes)
            Law being solved
        h: float
            The time step.
        qmf: QMF
            The filter whose parameters are
            being used in the simulation.
        ion: History
            The (sliced) array being used to
            determine the next time iteration's
            coordinates. It is not altered here,
            that is for the 'Corral' objects.

        Return
        ----------

        None
        """


        x = 1
        v = 2
        return x, v
    
class RK4(Generic_Solver):
    """
    This is the Runge-Kutta 4 Method
    """

    def __init__(self, specialClass=None):
        super().__init__(specialClass)
    
    def predict(self, law, h, qmf, ion):
        """
        Using the given law, integrate with the
        Runge-Kutta 4 method, to numerically
        solve these sets of coupled equations.

        K -> Coefficients for the Position Components (Depend on Velocity)
        G -> Coefficients for the Velocity Components (Depend on EoM)

        Parameters
        ----------

        law: Generic_Law (& sub-classes)
            Law being solved
        h: float
            The time step.
        qmf: QMF
            The filter whose parameters are
            being used in the simulation.
        ion: History
            The (sliced) array being used to
            determine the next time iteration's
            coordinates. It is not altered here,
            that is for the 'Corral' objects.

        Return
        ----------

        None
        """
        
        dZ = qmf.angFrequency * h

        
        # First Stage

        K1 = ion.vel
        G1 = law.compute(qmf, ion, ion.phaseTime[1])

        # Second Stage

        K2 = ion.vel + G1 * h / 2
        G2 = law.compute(qmf, ion + (K1 * h / 2, "X"), ion.phaseTime[1] + dZ / 2)

        # Third Stage

        K3 = ion.vel + G2 * h / 2
        G3 = law.compute(qmf, ion + (K2 * h / 2, "X"), ion.phaseTime[1] + dZ / 2)

        # Fourth Stage

        K4 = ion.vel + G3 * h
        G4 = law.compute(qmf, ion + (K3 * h, "X"), ion.phaseTime[1] + dZ)


        """
        End of the Runge-Kutta Method
        """
        
        x = ion.pos + h * (K1 + K4 + (K2 + K3) * 2) / 6
        v = ion.vel + h * (G1 + G4 + (G2 + G3) * 2) / 6  

        return x, v
    
class RK2(Generic_Solver):

    def __init__(self, specialClass=None):
        super().__init__(specialClass)
    
    def predict(self, law, h, qmf, ion):
        """
        Using the given law, integrate with the
        Runge-Kutta 4 method, to numerically
        solve these sets of coupled equations.

        K -> Coefficients for the Position Components (Depend on Velocity)
        G -> Coefficients for the Velocity Components (Depend on EoM)

        Parameters
        ----------

        law: Generic_Law (& sub-classes)
            Law being solved
        h: float
            The time step.
        qmf: QMF
            The filter whose parameters are
            being used in the simulation.
        ion: History
            The (sliced) array being used to
            determine the next time iteration's
            coordinates. It is not altered here,
            that is for the 'Corral' objects.

        Return
        ----------

        None
        """


        dZ = qmf.angFrequency * h
        
        # First Stage

        K1 = ion.vel
        G1 = law.compute(qmf, ion, ion.phaseTime[1])

        # Second Stage

        K2 = ion.vel + G1 * h / 2
        G2 = law.compute(qmf, ion + (K1 * h / 2, "X"), ion.phaseTime[1] + dZ / 2)
        
        x = ion.pos + h * K2
        v = ion.vel + h * G2

        return x, v
    
class Euler(Generic_Solver):
    """
    This is the Euler Method due to its simplicity and stability, 
    I will use it as a benchmark and for prototyping new versions
    of the Q521 program
    """

    def __init__(self, specialClass=None):
        super().__init__(specialClass)
    
    def predict(self, law, h, qmf, ion):
        """
        Using the given law, integrate with the
        Euler method, to numerically solve 
        these sets of coupled equations.

        Parameters
        ----------

        law: Generic_Law (& sub-classes)
            Law being solved
        h: float
            The time step.
        qmf: QMF
            The filter whose parameters are
            being used in the simulation.
        ion: History
            The (sliced) array being used to
            determine the next time iteration's
            coordinates. It is not altered here,
            that is for the 'Corral' objects.

        Return
        ----------

        None
        """
        
        v = ion.vel + h * law.compute(qmf, ion, ion.phaseTime[1])
        x = ion.pos + h * v

        return x, v