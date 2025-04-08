import numpy as np
import math
import sys

class Generic_Law:
    """
    This is a placeholder
    """

    def __init__(self, specialClass=None):
        self.special = specialClass
    
    def compute(self):
        # Also a placeholder, override this
        y = 0
        return y

class Mathieu(Generic_Law):
    """
    The Mathieu equations in the context of QMS physics.
    """
    def __init__(self, specialClass=None):
        super().__init__(specialClass)
    
    def compute(self, qmf, ion, phase):
        """
        Determine the acceleration caused by the electic field

        Parameters
        ----------

        qmf: QMF
            The filter whose parameters are
            being used in the simulation.
        ion: History
            The (sliced) array being used to
            determine the next time iteration's
            coordinates. It is not altered here,
            that is for the 'Corral' objects.
        phase: float
            The exact phase value being used in
            the wave form.

        Return
        ----------

        Numpy ndarray
            The acceleration for every ion at the
            time instant.
        """
        
        x = ion.pos
        v = ion.vel

        c0 = np.zeros([3,])
        c0[0] = -1
        c0[1] = 1

        c1 = ion.charge / (ion.mass * (qmf.inscRadius ** 2))
        c2 = qmf.directPotential + qmf.radioPotential * math.sin(phase)


        x = np.transpose(x)
        y = c0 * c1 * c2 * x
        y = np.transpose(y)
        
        return y

class Dawson_Entry(Generic_Law):
    """
    Fringing Field Model for a linear fields
    at the entry gap.
    """
    def __init__(self, specialClass=None):
        super().__init__(specialClass)
    
    def compute(self, qmf, ion, phase):
        """
        Determine the acceleration caused by a
        linear fringing field at the entry
        gap of the filter.

        Parameters
        ----------

        qmf: QMF
            The filter whose parameters are
            being used in the simulation.
        ion: History
            The (sliced) array being used to
            determine the next time iteration's
            coordinates. It is not altered here,
            that is for the 'Corral' objects.
        phase: float
            The exact phase value being used in
            the wave form.

        Return
        ----------

        Numpy ndarray
            The acceleration for every ion at the
            time instant.
        """
        
        x = ion.pos
        v = ion.vel
        gap = qmf.entryGap

        c0 = np.zeros([3,])
        c0[0] = -1
        c0[1] = 1

        c1 = ion.charge / (ion.mass * (qmf.inscRadius ** 2))
        c2 = qmf.directPotential + qmf.radioPotential * math.sin(phase)
        c3 = np.zeros([3, ion.number])
        c3[0:2, :] = (x[2, :]) / gap

        x = np.transpose(x)
        c3 = np.transpose(c3)
        y = c0 * c1 * c2 * (c3 * x)
        y = np.transpose(y)
        
        return y 

class Dawson_Exit(Generic_Law):
    """
    Fringing Field Model for a linear fields
    at the exit gap.
    """

    def __init__(self, specialClass=None):
        super().__init__(specialClass)
    
    def compute(self, qmf, ion, phase):
        """
        Determine the acceleration caused by a
        linear fringing field at the exit
        gap of the filter.

        Parameters
        ----------

        qmf: QMF
            The filter whose parameters are
            being used in the simulation.
        ion: History
            The (sliced) array being used to
            determine the next time iteration's
            coordinates. It is not altered here,
            that is for the 'Corral' objects.
        phase: float
            The exact phase value being used in
            the wave form.

        Return
        ----------

        Numpy ndarray
            The acceleration for every ion at the
            time instant.
        """
        
        x = ion.pos
        v = ion.vel
        gap = qmf.exitGap
        start = qmf.entryGap + qmf.span

        c0 = np.zeros([3,])
        c0[0] = -1
        c0[1] = 1

        c1 = ion.charge / (ion.mass * (qmf.inscRadius ** 2))
        c2 = qmf.directPotential + qmf.radioPotential * math.sin(phase)
        c3 = np.zeros([3, ion.number])
        c3[0:1, :] = 1 - ((x[2, :] - start) / gap)


        x = np.transpose(x)
        c3 = np.transpose(c3)
        y = c0 * c1 * c2 * (c3 * x)
        y = np.transpose(y)
        
        return y

class HM_Entry(Generic_Law):
    """
    Hunter-McIntosh model of fringing fields
    at the entry gap.
    """

    def __init__(self, specialClass=None):
        super().__init__(specialClass)
    
    def compute(self, qmf, ion, phase):
        """
        Determine the acceleration caused by a
        Hunter-McIntosh fringing field 
        at the entry gap of the filter.

        Parameters
        ----------

        qmf: QMF
            The filter whose parameters are
            being used in the simulation.
        ion: History
            The (sliced) array being used to
            determine the next time iteration's
            coordinates. It is not altered here,
            that is for the 'Corral' objects.
        phase: float
            The exact phase value being used in
            the wave form.

        Return
        ----------

        Numpy ndarray
            The acceleration for every ion at the
            time instant.
        """
        
        x = ion.pos
        v = ion.vel
        gap = qmf.entryGap

        c0 = np.zeros([3,])
        c0[0] = -1
        c0[1] = 1

        c1 = ion.charge / (ion.mass * (qmf.inscRadius ** 2))
        c2 = qmf.directPotential + qmf.radioPotential * math.sin(phase)
        c3 = np.zeros([3, ion.number])
        z = (x[2, :]) / gap
        c3[0:2, :] = 1 - np.exp( (0- qmf.a1 * z - qmf.b1 * (z ** 2)))

        # Beta-only
        #save(c3[1, 0])
        # Close

        x = np.transpose(x)
        c3 = np.transpose(c3)
        y = c0 * c1 * c2 * (c3 * x)
        y = np.transpose(y)
        
        return y

class HM_Exit(Generic_Law):
    """
    (Mirrored) Hunter-McIntosh model of fringing fields
    at the exit gap.
    """

    def __init__(self, specialClass=None):
        super().__init__(specialClass)
    
    def compute(self, qmf, ion, phase):
        """
        Determine the acceleration caused by a
        Hunter-McIntosh fringing field 
        at the exit gap of the filter.

        Parameters
        ----------

        qmf: QMF
            The filter whose parameters are
            being used in the simulation.
        ion: History
            The (sliced) array being used to
            determine the next time iteration's
            coordinates. It is not altered here,
            that is for the 'Corral' objects.
        phase: float
            The exact phase value being used in
            the wave form.

        Return
        ----------

        Numpy ndarray
            The acceleration for every ion at the
            time instant.
        """
        
        x = ion.pos
        v = ion.vel
        gap = qmf.entryGap

        c0 = np.zeros([3,])
        c0[0] = -1
        c0[1] = 1

        c1 = ion.charge / (ion.mass * (qmf.inscRadius ** 2))
        c2 = qmf.directPotential + qmf.radioPotential * math.sin(phase)
        c3 = np.zeros([3, ion.number])
        z = (x[2, :]) / gap
        c3[0:2, :] = 1 - np.exp( (0- qmf.a1 * z - qmf.b1 * (z ** 2)))

        # Beta-only
        #save(c3[1, 0])
        # Close

        x = np.transpose(x)
        c3 = np.transpose(c3)
        y = c0 * c1 * c2 * (c3 * x)
        y = np.transpose(y)
        
        return y


"""
Monitoring Functions
"""

def save(number):
    with open("amp.txt", "a") as file:
        file.write(str(number) + "\n")