import os
import math
import time
from datetime import datetime
from copy import deepcopy
import numpy as np
from matplotlib import pyplot as plt
from qmf import QMF
from history import History, Cluster, History
from corral import Corral
from law import Mathieu, Dawson_Entry, Dawson_Exit, HM_Entry
from solver import Euler, RK4

"""
Simulator Class Definition
"""

class Simulator:
    """
    Q521's simulation routines wrapped in a single object type.
    """
    def __init__(self):
        """
        Initialise the various attributes of the Simulator object.

        Parameters
        ----------

        None

        Return
        ----------

        None
        """
        
        # Custom Objects
        self.qmf = None
        self.ion = None
        self.skimIon = None
        self.queue = None

        # Standard Objects
        self.h = None
        self.div = 100 # Default -> 1/ 100th of a RF Period as the duration of the phase-time step

        self.qmfString = ""
        self.ionString = ""

        self.qmfTag = ""
        self.ionTag = ""

        self.lap1 = None
        self.lap2 = None
        self.lap3 = None
        
        # Display-oriented Objects
        self.fig1 = None
        self.fig2 = None
        self.fig3 = None
        self.fig4 = None

    def load(self):
        """
        Creates the 'History' and 'QMF' object instances and loads in their
        associated parameters from the ion and filter directories.
        Note that before this method is called, a file name and tag
        must be assigned to the 'Simulator' object's qmfString,
        qmfTag, ionString and ionTag respectively.

        Parameters
        ----------

        None

        Return
        ----------

        None
        """
        # Create instances based on the chosen file and tag
        self.qmf = QMF(self.qmfString, self.qmfTag)
        self.ion = Cluster(self.ionString)

        # Load the data into the object
        self.qmf.load_csv()
        self.ion.load_csv()
    
    def setup(self):
        """
        Determines the number of steps for each quadrupole
        region as well as pre-allocating the memory space
        for the Numpy array attributes of the 'History'
        object. The 'Corral' objects are also assigned to
        the 'queue' attribute (a list of usually 3 objects).
        IT IS ASSUMED THAT THE 'load' METHOD HAS ALREADY BEEN
        CALLED.

        Parameters
        ----------

        None

        Return
        ----------

        None
        """
        
        # ASSUME that load() has been called
        
        # Time Step Definition
        self.h = 1 / (self.div * self.qmf.pureFrequency)

        # Step Lengths -> Number of Steps (Needed for the Corral Queue)
        self.lap1 = math.ceil((self.qmf.entryGap / self.qmf.injSpeed) / self.h)
        self.lap2 = math.ceil((self.qmf.span / self.qmf.injSpeed) / self.h)
        self.lap3 = math.ceil((self.qmf.exitGap / self.qmf.injSpeed) / self.h)

        print(self.lap1, self.lap2, self.lap3)

        # Ion Setup -> Note that the number of elements is the sum of 
        # the num of steps + 1 for the initial
        self.ion.setup(1 + self.lap1 + self.lap2 + self.lap3, self.qmf)

        if isinstance(self.ion, History):
            pass
        elif isinstance(self.ion, Cluster):
            G = len(self.ion)
            x, y = (0, 0)
            for g in range(G):
                _, t, n = self.ion[g].pos.shape
                x += n
                y = t
            print(f"The dimensions of the total array are: {y} time points and {x} ions")

            print( isinstance(self.ion, Cluster) )

        # Queue Creation
        self.queue = [
            Corral(Euler, HM_Entry, self.h, 0, self.lap1, self.qmf),
            Corral(Euler, Mathieu, self.h, self.lap1, self.lap2, self.qmf),
            Corral(Euler, Dawson_Exit, self.h, self.lap2, self.lap3, self.qmf),
            ]

    def run(self):
        """
        Performs the simulation across all the elements
        in the 'Corral' queue.

        Parameters
        ----------

        None

        Return
        ----------

        None
        """

        # Loop through each operation
        for corral in self.queue:
            corral.lockOn(self.ion)
        
        # Keep track of current time
        pureDate = datetime.now()
        s = pureDate.strftime("%Y-%m-%d %H_%M_%S")
        self.date = s[0:10]
        self.hour = s[11:19]

        print("Simulation is Complete!")
    
    def clean(self):
        """
        Searches the 'ion' attribute (a 'History' object)
        for any lost ions and eliminates them from a duplicate
        of the attribute called 'skimIon'. The 'skimIon' 
        attribute is later used for creating plots without any
        divergent curves.

        Parameters
        ----------

        None

        Return
        ----------

        None
        """

        if isinstance(self.ion, History):
            
            # Loop through each ion for each time step ---> {Not very efficient, but should be fine}
            # If ion is lost, mark its index as false
            # in the 'status' array
            _, T, N = self.ion.pos.shape
            for t in range(T):
                for n in range(N):
                    C = self.ion.status[n]
                    if not C:
                        pass
                    else:
                        # Assign this iteration's state arrays
                        x, y = self.ion.pos[0:2, t, n]
                        r = self.qmf.inscRadius
                        
                        # Determine the difference between u and r_o
                        d1 = abs(x) - r
                        d2 = abs(y) - r
                        
                        # Tick as false or pass if lost
                        if (d1 >=0) or (d2 >= 0):
                            self.ion.status[n] = 0
                        else:
                            pass
            
            # Create a list with all the indices of the lost ions
            lst = np.where(self.ion.status == 0)[0]
            lst = np.flip(lst) # flip it, so that we can remove the laters indices first
            
            # Determine the transmission rate
            p = 100 * (N - lst.shape[0]) / N
            print([lst, lst.shape, f"Transmission Rate: {p}%"])

            # Create the 'skimIon' attribute from a replicate of 'ion'
            self.skimIon = self.ion.replicate()

            # Delete the "card" array because its lost
            for i in lst:
                self.skimIon.pos = np.delete(self.skimIon.pos, i, 2)
                self.skimIon.vel = np.delete(self.skimIon.vel, i, 2)
        
        elif isinstance(self.ion, Cluster):
            
            self.skimIon  = deepcopy(self.ion)
            G = len(self.skimIon)

            for g in range(G):
                
                _, T, N = self.skimIon._values[g].pos.shape
                
                for t in range(T):
                    
                    for n in range(N):
                        
                        C = self.skimIon._values[g].status[n]
                        
                        if not C:
                            pass
                        
                        else:
                            # Assign this iteration's state arrays
                            x, y = self.skimIon._values[g].pos[0:2, t, n]
                            r = self.qmf.inscRadius
                            
                            # Determine the difference between u and r_o
                            d1 = abs(x) - r
                            d2 = abs(y) - r
                            
                            # Tick as false or pass if lost
                            if (d1 >= 0) or (d2 >= 0):
                                self.skimIon._values[g].status[n] = 0
                            else:
                                pass
            
            # loop through all history objects in the container
            # remove the lost ions
            for g in range(G):
                
                A = self.skimIon._values[g].status == 0
                B = np.where(A)
                tag = self.skimIon._values[g].tag
                nP = self.skimIon._values[g].number

                lst = B[0] # possible bug point
                lst = np.flip(lst) # flip it
                
                # Determine the transmission rate
                p = 100 * (nP - lst.shape[0]) / nP
                a = self.skimIon._values[g].a
                q = self.skimIon._values[g].q
                print([lst, lst.shape, f"Transmission Rate for this Ion Species {tag}: {p}%"])


                # Delete the "card" array because its lost
                for i in lst:
                    self.skimIon._values[g].pos = np.delete(self.skimIon._values[g].pos, i, 2)
                    self.skimIon._values[g].vel = np.delete(self.skimIon._values[g].vel, i, 2)
                _, _, c = self.skimIon._values[g].pos.shape
                print(f"Number of ions still left: {c} for species {g}")
    
    def mass_scan(self):
        if isinstance(self.ion, History):
            pass
        elif isinstance(self.ion, Cluster):
            I = 1
            G = len(self.skimIon)
            v = np.ndarray([G,])
            u = np.zeros([G,])
            

            # This loop creates the horizontal array for the spectrum
            for g in range(G):
                v[g] = (self.skimIon._values[g].mass * 1.66e+27) / abs(self.skimIon._values[g].charge * 1.6e+19)
            
            # This nested loop iterates over mass scans and adds the number
            # of transmitted ions to the total count
            for i in range(I):
                for g in range(G):
                    _, _, n = self.skimIon._values[g].pos.shape
                    u[g] += n

            plt.bar(v, u, width=0.025)
            plt.show()


    def temp_display(self):
        """
        {TEMPORARY} Creates the trajectory curves 
        as well as the phase curves of the multiple
        ions in the 'skimIon' attribute.

        Parameters
        ----------

        None

        Return
        ----------

        None
        """
        try:
            if isinstance(self.ion, History):
                self.hstTime = self.skimIon.phaseTime[0, :]
                self.hstTrack = self.skimIon.pos[2, :, 0]

                # Create the Matplotlib figures
                self.fig1, (self.ax1_1, self.ax1_2) = plt.subplots(2, 1)
                self.fig2, (self.ax2_1, self.ax2_2) = plt.subplots(2, 1)
                self.fig3, (self.ax3_1, self.ax3_2) = plt.subplots(2, 1)
                self.fig4, (self.ax4_1, self.ax4_2) = plt.subplots(2, 1)

                # Plot the trajectory curves in Figure 1
                for i in range(self.skimIon.number):
                    self.ax1_1.plot(self.hstTrack, self.skimIon.pos[0, :, i])
                    self.ax1_2.plot(self.hstTrack, self.skimIon.pos[1, :, i])

                # Plot the trajectory curves in Figure 2
                for i in range(self.skimIon.number):
                    self.ax2_1.plot(self.skimIon.pos[0, :, i], self.skimIon.vel[0, :, i] / 500)
                    self.ax2_2.plot(self.skimIon.pos[1, :, i], self.skimIon.vel[1, :, i] / 500)   

                # Plot the trajectory curves in Figure 3
                x = self.skimIon.pos[0, 0, :]
                y = self.skimIon.pos[1, 0, :]

                x_f = self.skimIon.pos[0, -1, :]
                y_f = self.skimIon.pos[1, -1, :]
                
                for i in range(self.skimIon.number):
                    
                    self.ax3_1.scatter(x, y)
                    self.ax3_2.scatter(x_f, y_f)

                # Plot the trajectory curves in Figure 4
                shortTrack = self.skimIon.pos[2, 0 : self.lap1 + 1, 0]
                
                for i in range(self.skimIon.number):
                    self.ax4_1.plot(shortTrack, self.skimIon.pos[0, 0 : self.lap1 + 1, i])
                    self.ax4_2.plot(shortTrack, self.skimIon.pos[1, 0 : self.lap1 + 1, i])
                
                # Axis Aesthetics
                self.ax1_1.set_title("XZ Space")
                self.ax1_2.set_title("YZ Space")

                self.ax2_1.set_title("XV Space")
                self.ax2_2.set_title("YV Space")

                self.ax3_1.set_title("Initial Position")
                self.ax3_2.set_title("Final Position")

                self.ax4_1.set_title("Fringe XZ Space")
                self.ax4_2.set_title("Fringe YZ Space")

                # Create results directory, if missing
                os.makedirs("results", exist_ok=True)
                os.makedirs("results/" + self.date, exist_ok=True)

                # Store the figures 
                self.ion.save("results/" + self.date, self.date, self.hour)

                self.fig1.savefig("results/"+ self.date + "/SpaceCurve " + self.date + " " + self.hour + ".svg", transparent=True)
                self.fig2.savefig("results/"+ self.date + "/PhaseCurve " + self.date + " " + self.hour + ".svg", transparent=True)
                self.fig3.savefig("results/"+ self.date + "/HeatMap " + self.date + " " + self.hour + ".svg", transparent=True)
                self.fig4.savefig("results/"+ self.date + "/FringeCurve " + self.date + " " + self.hour + ".svg", transparent=True)
                plt.tight_layout()
                plt.show()

            elif isinstance(self.ion, Cluster):
                G = len(self.ion)

                self.hstTime = self.ion._values[0].phaseTime[0, :]
                self.hstTrack = self.ion._values[0].pos[2, :, 0]

                # Create the Matplotlib figures
                self.fig1, (self.ax1_1, self.ax1_2) = plt.subplots(2, 1)
                self.fig2, (self.ax2_1, self.ax2_2) = plt.subplots(2, 1)
                self.fig3, (self.ax3_1, self.ax3_2) = plt.subplots(2, 1)
                self.fig4, (self.ax4_1, self.ax4_2) = plt.subplots(2, 1)


                # Plot the trajectory curves in Figure 1
                for g in range(G):
                    _, _, N = self.skimIon._values[g].pos.shape
                    for i in range(N):
                        self.ax1_1.plot(self.hstTrack, self.skimIon._values[g].pos[0, :, i])
                        self.ax1_2.plot(self.hstTrack, self.skimIon._values[g].pos[1, :, i])


                # Plot the trajectory curves in Figure 2
                for g in range(G):
                    _, _, N = self.skimIon._values[g].pos.shape
                    for i in range(N):
                        self.ax2_1.plot(self.skimIon._values[g].pos[0, :, i], self.skimIon._values[g].vel[0, :, i] / 500)
                        self.ax2_2.plot(self.skimIon._values[g].pos[1, :, i], self.skimIon._values[g].vel[1, :, i] / 500)

                print("Past")

                # Plot the trajectory curves in Figure 3
                for g in range(G):
                    
                    x = self.skimIon._values[g].pos[0, 0, :]
                    y = self.skimIon._values[g].pos[1, 0, :]

                    x_f = self.skimIon._values[g].pos[0, -1, :]
                    y_f = self.skimIon._values[g].pos[1, -1, :]
                        
                    self.ax3_1.scatter(x, y)
                    self.ax3_2.scatter(x_f, y_f)


                print("Past Perfect")



                # Plot the trajectory curves in Figure 4
                shortTrack = self.ion._values[0].pos[2, 0 : self.lap1 + 1, 0]
                
                for g in range(G):
                    _, _, N = self.skimIon._values[g].pos.shape
                    for i in range(N):
                        
                        self.ax4_1.plot(shortTrack, self.skimIon._values[g].pos[0, 0 : self.lap1 + 1, i])
                        self.ax4_2.plot(shortTrack, self.skimIon._values[g].pos[1, 0 : self.lap1 + 1, i])

                print("Past More-than-perfect")
                
                # Axis Aesthetics
                self.ax1_1.set_title("XZ Space")
                self.ax1_2.set_title("YZ Space")

                self.ax2_1.set_title("XV Space")
                self.ax2_2.set_title("YV Space")

                self.ax3_1.set_title("Initial Position")
                self.ax3_2.set_title("Final Position")

                self.ax4_1.set_title("Fringe XZ Space")
                self.ax4_2.set_title("Fringe YZ Space")

                # Create results directory, if missing
                os.makedirs("results", exist_ok=True)
                os.makedirs("results/" + self.date, exist_ok=True)

                # Store the figures 
                self.ion.save("results/" + self.date, self.date, self.hour)

                self.fig1.savefig("results/"+ self.date + "/SpaceCurve " + self.date + " " + self.hour + ".svg", transparent=True)
                self.fig2.savefig("results/"+ self.date + "/PhaseCurve " + self.date + " " + self.hour + ".svg", transparent=True)
                self.fig3.savefig("results/"+ self.date + "/HeatMap " + self.date + " " + self.hour + ".svg", transparent=True)
                self.fig4.savefig("results/"+ self.date + "/FringeCurve " + self.date + " " + self.hour + ".svg", transparent=True)
                plt.tight_layout()
                plt.show()
        except:
            pass

sim = Simulator()
sim.qmfString = "filter1.csv"
sim.qmfTag = "McIntosh Filter"
sim.ionString = "ionspec2.csv"
sim.ionTag = "Sodium"


startTime = time.time()

sim.load()
sim.setup()
sim.run()

endTime = time.time()

print(f"The simulation took: {endTime - startTime} seconds")

sim.clean()
sim.temp_display()
sim.mass_scan()
print("End of the Line")
