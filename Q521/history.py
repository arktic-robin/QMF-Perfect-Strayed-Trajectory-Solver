import numpy as np
import csv
from copy import deepcopy
import random as rng

class History:
    """
    This class is intended to store the parameters and coordinates
    of a specific grouping of ions through time or phase. There are
    two main ways to declare this object type:
    1) Specify the (file)name and tag and then use the load_csv method
    2) Specify a filename, tag and all the parameters [Less Accesible]
    """
    
    def __init__(self, name, tag, param=[]):
        """
        Initialise the filter object with a csv file
        or through the manual definition of the QMF
        parameters.

        Parameters
        ----------

        name: string
            The csv file to read the parameters
            from.
        tag: string
            The filter to be used from the QMF
            set in the csv file.
        param: list
            {NOT OPTIMAL} Manual definition of
            the parameters of the quadrupole.
            The name and tag must be set to
            None.

        Return
        ----------

        None
        """
        
        # Input Parameters
        self.name = name
        self.tag = tag
        
        # Manual Definition
        if name == None:
            self.tag = param[0]
            self.number = param[1]
            self.mass = param[2]
            self.charge = param[3]
            self.spX = param[4]
            self.spY = param[5]
            self.spVX = param[6]
            self.spVY = param[7]

    def setup(self, nodes, filter):
        """
        Pre-allocate the time and coordinate 
        arrays and creates the initial conditions
        of the simulation for every single ion.

        Parameters
        ----------

        nodes: int
            Number of iterations which
            the simulation will occupy.
        filter: QMF
            The filter which will be used
            for the simulation.

        Return
        ----------

        None
        """
        
        self.phaseTime = np.zeros([2, nodes])
        self.status = np.ones([self.number,])
        self.pos = np.zeros([3, nodes, self.number])
        self.vel = np.zeros([3, nodes, self.number])

        self.phaseTime[1, 0] = filter.iniPhase
        for i in range(self.number):
            del_x = self.spX * (2 * rng.random() - 1)
            del_y = self.spY * (2 * rng.random() - 1)

            del_vx = self.spVX * (2 * rng.random() - 1)
            del_vy = self.spVY * (2 * rng.random() - 1)

            self.pos[0, 0, i] = del_x
            self.pos[1, 0, i] = del_y
            self.pos[2, 0, i] = 0

            self.vel[0, 0, i] = del_vx
            self.vel[1, 0, i] = del_vy
            self.vel[2, 0, i] = filter.injSpeed
    
    def burn(self):
        """
        Pre-allocates a trimmed array intended
        for a single time step.

        Parameters
        ----------

        None

        Return
        ----------

        None
        """
        
        self.phaseTime = np.zeros([2, 1])
        self.pos = np.zeros([3, self.number])
        self.vel = np.zeros([3, self.number])
    
    def load_csv(self):
        """
        Read the ion species set 
        and reassign the values
        of the ion attributes.

        Parameters
        ----------

        None

        Return
        ----------

        None
        """

        with open("ion/" + self.name) as f:
            rows = csv.reader(f, delimiter=",")
            for row in rows:
                if row[0] == self.tag:
                    self.number = int(row[1])
                    self.mass = float(row[2]) * 1.66e-27
                    self.charge = float(row[3]) * 1.6e-19
                    self.spX = float(row[4]) * 1e-3
                    self.spY = float(row[5]) * 1e-3
                    self.spVX = float(row[6]) * 1e-3
                    self.spVY = float(row[7]) * 1e-3

    def eject(self, j):
        """
        {Legacy} THIS IS VERY INEFFICIENT, FIND A BETTER SOLUTION
        """
        clone = self.replicate()
        clone.phaseTime = clone.phaseTime[:, j]
        clone.pos = clone.pos[:, j, :]
        clone.vel = clone.vel[:, j, :]
        return clone

    def fastEject(self, j):
        """
        Splits of a slice of the array at time
        index J, unlinking the arrays in the
        memory of the computer.

        Parameters
        ----------

        j: int
            Time index.

        Return
        ----------

        History
            A trimmed 'History' object with
            only the coordinates at time J.
        """
        
        # Duplicates a new instance of the History object
        y =  History(None, None, param=[
            self.tag,
            self.number,
            self.mass,
            self.charge,
            self.spX,
            self.spY,
            self.spVX,
            self.spVY
        ])

        # Creates an matrix slice attribute and inserts Current value
        y.phaseTime = np.zeros([2,])
        y.pos = np.zeros([3, y.number])
        y.vel = np.zeros([3, y.number])
        
        
        y.phaseTime = y.phaseTime + self.phaseTime[:, j]
        y.pos = y.pos + self.pos[:, j, :]
        y.vel = y.vel + self.vel[:, j, :]

        return y
    
    def replicate(self):
        """
        {Legacy} THIS IS VERY INEFFICIENT, FIND A BETTER SOLUTION
        """
        
        clone = deepcopy(self)
        return clone
    
    def save(self, dir, date, hour):
        """
        Save the coordinate arrays in a txt/csv file.
        The structure of the file is, N * [x,y,z] for
        columns and the number of rows is the total
        number of time points.

        Parameters
        ----------

        None

        Return
        ----------

        None
        """

        with open(dir + "/set "+ date + " " + hour + ".csv", "w") as file:
            _, I, J = self.pos.shape
            
            for i in range(I):
                line = ""
                for j in range(J):
                    line = line + str(self.pos[0, i, j]) + ","
                    line = line + str(self.pos[1, i, j]) + ","
                    line = line + str(self.pos[2, i, j]) + ","
                file.write(line+ "\n")

    
    def step(self, j, dt, w):
        """
        Determine the next time point's phase
        and time values.

        Parameters
        ----------

        j: int
            Time index.
        dt: float
            Time step
        w: float
            Angular frequency.

        Return
        ----------

        None
        """
        
        self.phaseTime[0, j] = self.phaseTime[0, j - 1] + dt
        self.phaseTime[1, j] = self.phaseTime[1, j - 1] + w * dt
    
    def __add__(self, other):
        """
        Addition Operation
        """
        # Certain Things
        x = other[0]
        branch = other[1]


        # Duplicates a new instance of the History object
        y =  History(None, None, param=[
            self.tag,
            self.number,
            self.mass,
            self.charge,
            self.spX,
            self.spY,
            self.spVX,
            self.spVY
        ])

        # Creates an matrix slice attribute and inserts Current value
        y.burn()
        y.pos = y.pos + self.pos
        y.vel = y.vel + self.vel
        

        match branch:
            case "X":
                y.pos = y.pos + x
            case "V":
                y.vel = y.vel + x
        return y

class Cluster:
    """
    This class is intended to store multiple 'History'
    objects, all being created and set-up from a
    pre-defined parameter set (csv).
    """
    
    def __init__(self, name):        
        # Input Parameters
        self.name = name
        self._values = [] # Private attribute
    
    def load_csv(self):
        
        with open("ion/" + self.name) as f:
            
            i = -2
            rows = csv.reader(f, delimiter=",")
            

            for row in rows:
                i += 1
                if row[0] == "Tag":
                    continue

                row[1] = int(row[1])
                row[2] = float(row[2]) * 1.66e-27
                row[3] = float(row[3]) * 1.6e-19
                row[4] = float(row[4]) * 1e-3
                row[5] = float(row[5]) * 1e-3
                row[6] = float(row[6]) * 1e-3
                row[7] = float(row[7]) * 1e-3
                
                ion = History(None, None, param=row)
                self.append(ion)
                print(row)
    
    def setup(self, nodes, filter):
        """
        Setup for the Cluster Type
        """
        G = len(self)

        for g in range(G):
        
            den = self._values[g].mass * (filter.angFrequency * filter.inscRadius) ** 2
            self._values[g].a = (4 * self._values[g].mass * filter.directPotential) / den
            self._values[g].q = 2 * self._values[g].mass * filter.radioPotential / den
            
            self._values[g].phaseTime = np.zeros([2, nodes])
            self._values[g].status = np.ones([self._values[g].number,])
            self._values[g].pos = np.zeros([3, nodes, self._values[g].number])
            self._values[g].vel = np.zeros([3, nodes, self._values[g].number])

            self._values[g].phaseTime[1, 0] = filter.iniPhase
            for i in range(self._values[g].number):
                del_x = self._values[g].spX * (2 * rng.random() - 1)
                del_y = self._values[g].spY * (2 * rng.random() - 1)

                del_vx = self._values[g].spVX * (2 * rng.random() - 1)
                del_vy = self._values[g].spVY * (2 * rng.random() - 1)

                self._values[g].pos[0, 0, i] = del_x
                self._values[g].pos[1, 0, i] = del_y
                self._values[g].pos[2, 0, i] = 0

                self._values[g].vel[0, 0, i] = del_vx
                self._values[g].vel[1, 0, i] = del_vy
                self._values[g].vel[2, 0, i] = filter.injSpeed
    
    def append(self, val):
        self._values.append(val)

    def save(self, dir, date, hour):
        pass
  
    def step(self, j, dt, w):
        pass
    
    def __getitem__(self, i):
        return self._values[i]
    
    def __setitem__(self, i, val):
        self._values[i] = val
    
    def __len__(self):
        return len(self._values)
