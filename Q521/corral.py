from history import History, Cluster

"""
Corral Object
"""

class Corral:
    """
    Writes onto a 'History' object the coordinates
    and velocities of ions through specific field
    conditions (like inside the span, or fringes).
    """

    def __init__(self, solver, law, timeStep, startIndex, steps, filter):
        """
        Initialises the 'Corral' object

        Parameters
        ----------

        solver: Generic_Solver (& sub-classes)
            The type of solver to be used in the
            corral's calculations (Euler/RK4/RK2).
            Note that the class is the input, not
            an instance.
        
        law: Generic_Law (& sub-classes)
            The law being solved, be it Mathieu,
            Dawson or Hunter-McIntosh (HM).
            Note that the class is the input, not
            an instance; Also, for fringing field
            models, specify whether the law is to
            be used for the entry gap or exit,
            that is, 'HM_Entry' or 'Dawson_Exit'.
        
        timeStep: float
            The time step used for teh simulation,
            ideally use a time step of T/k, with 
            k >=100, and T being the radio-period.
        
        startIndex: int
            The time index to start writing from
            in the 'History' object being written
            on.
        
        steps: int
            Quantity of time steps to be run.

        filter: QMF
            The filter which is being used.

        Return
        ----------

        None
        """
        
        # Object Assignments
        self.h = timeStep
        self.start  = startIndex
        self.steps = steps 
        self.qmf = filter

        # Object Constructor Calls
        self.solver = solver()
        self.law = law()
    
    def lockOn(self, ion):
        """
        Start running calculations and writing
        the coordinate data onto the input 'ion'
        object.

        Parameters
        ----------

        ion: History or Cluster
            The object which is being written to.
            For 'Cluster' objects, the program will
            perform the simulation for each individual
            'History' object in the container.

        Return
        ----------

        None
        """

        # Determine whether this object isn't operable
        # else, run the calculation for each time step
        if self.steps == None or self.steps == 0:
            pass
        
        else:
            # Create iterable list of each time step (index - 1)
            I = list(range(self.start, self.start + self.steps))

            # If used in "Dive Mode", run the calculation as normal
            if isinstance(ion, History):

                for i in I:
                    i += 1
                    #print(f"Current Stage of the Corraling Process --------------> {100 * (i - I[0] - 1) / (I[-1] - I[0])}%")

                    # Perform the numerical integration  
                    ion.step(i, self.h, self.qmf.angFrequency) # Determine next index's time value
                    tempIon = ion.fastEject(i - 1) # Extract a i-time slice of the array
                    x, v = self.solver.predict(self.law, self.h, self.qmf, tempIon) # Use the solver to determine future values

                    # Write these to the pre-allocated coordinate slots
                    ion.pos[:, i, :] += x
                    ion.vel[:, i, :] += v
            
            # If used in "Scan Mode", run the calculation for each
            # 'History' object
            elif isinstance(ion, Cluster):
                # Get the Length of the 'Cluster' object ie Number of Ion Species
                G = len(ion)
                G = list(range(G))

                # Iterate over each ion species
                for g in G:
                    
                    # Notice that the progress meter is based on the number of completed Ion Species now
                    try:
                        print(f"Current Stage of the Cluster Corraling Process --------------> {100 * (g - G[0]) / (G[-1] - G[0])}%")
                    except:
                        pass
                    
                    # Iterate over each time point
                    for i in I:
                        i += 1
                        
                        # Perform the numerical integration for the 
                        ion[g].step(i, self.h, self.qmf.angFrequency) # Note that we are referencing  cluster[g]
                        tempIon = ion[g].fastEject(i - 1)
                        x, v = self.solver.predict(self.law, self.h, self.qmf, tempIon)

                        # Write these to the pre-allocated coordinate slots
                        # of the appropriate g-th 'History' object
                        ion[g].pos[:, i, :] += x
                        ion[g].vel[:, i, :] += v
