import numpy as np
from numpy import ones, zeros
import random as rng
import csv
from datetime import datetime
from time import sleep
import os
from matplotlib import pyplot as plt

__all__ = ["State"]

class State:
    
    def __init__(self):

        self.tag = None
        
        self.time = None
        self.state = None
        self.history = None
        self.solver = None
        self.h = None
        
        self.length = None
        self.rad = None

        self.number = None
        self.mass = None
        self.charge = None
        self.sprX = None
        self.sprY = None
        self.sprVX = None
        self.sprVY = None
        self.sprVZ = None
        self.injSpeed = None
        self.iniPhase = None
    
    def build(self, timepoints):

        self.time = zeros([timepoints])
        self.state = zeros([7, self.number])
        self.history = zeros([7, self.number, timepoints])

        # Determine time values
        for t in range(timepoints - 1):
            tempo = t * self.h
            self.time[t] = tempo

        # Set ion initial conditions
        for n in range(self.number):
            
            delX = self.sprX * (2 * rng.random() - 1)
            delY = self.sprY * (2 * rng.random() - 1)
            
            delVX = self.sprVX * (2 * rng.random() - 1)
            delVY = self.sprVY * (2 * rng.random() - 1)
            delVZ = self.sprVZ * rng.random() * 0.001

            self.state[0, n] = 1
            self.state[1:4, n] += [delX, delY, 0]
            self.state[4:7, n] += [delVX, delVY, self.injSpeed + delVZ]
    
    def load(self, name):
        
        with open(name) as f:
            
            rows = csv.reader(f, delimiter=",")

            for row in rows:
                
                if row[0] == self.tag:
                    
                    self.number = int(row[1])
                    self.mass = float(row[2]) * 1.66e-27
                    self.charge = float(row[3]) * 1.6e-19
                    self.sprX = float(row[4]) * 1e-3
                    self.sprY = float(row[5]) * 1e-3
                    self.sprVX = float(row[6]) * 1e-3
                    self.sprVY = float(row[7]) * 1e-3
                    self.sprVZ = float(row[8]) * 1e3
                    self.injSpeed = float(row[9]) * 1e3
                    self.iniPhase = float(row[10])
    
    def clock(self):
        # Keep track of current time
        pureDate = datetime.now()
        s = pureDate.strftime("%Y-%m-%d %H_%M_%S")
        self.date = s[0:10]
        self.hour = s[11:19]

        return self.date, self.hour
    
    def save(self, groupCode=""):

        (date, time) = self.clock()
        # Create results directory, if missing
        os.makedirs("results", exist_ok=True)
        os.makedirs("results/" + date, exist_ok=True)
        
        with open( "results/"+ date + "/" + "ion_history " + groupCode + " " + time + ".csv", "w") as file:

            _, N, T = self.history.shape

            # Header loop

            line = ""
            for n in range(N):
                line = line + "L" + str(n) + ","
                line = line + "X" + str(n) + ","
                line = line + "Y" + str(n) + ","
                line = line + "Z" + str(n) + ","
                line = line + "VX" + str(n) + ","
                line = line + "VY" + str(n) + ","
                line = line + "VZ" + str(n) + ","
            
            line = line[0:-2]

            file.write(line+"\n")

            
            # Data loop
            for t in range(T):

                line = ""
                
                for n in range(N):
                    line = line + str(self.history[0, n, t]) + ","
                    line = line + str(self.history[1, n, t]) + ","
                    line = line + str(self.history[2, n, t]) + ","
                    line = line + str(self.history[3, n, t]) + ","
                    line = line + str(self.history[4, n, t]) + ","
                    line = line + str(self.history[5, n, t]) + ","
                    line = line + str(self.history[6, n, t]) + ","
                
                line = line[0:-2]
                
                file.write(line+"\n")

    def display(self, save_fig=True, groupCode=""):

        l = self.length
        r = self.rad

        N = self.number

        t = self.time

        x = None
        y = None
        z = None

        u = None
        v = None
        w = None


        # Create the Matplotlib figures
        self.fig1, (self.ax1_1, self.ax1_2) = plt.subplots(2, 1)
        self.fig2, (self.ax2_1, self.ax2_2, self.ax2_3) = plt.subplots(3, 1)
        self.fig3, (self.ax3_1, self.ax3_2) = plt.subplots(2, 1)

        for n in range(N):

            # Assignments

            x = self.history[1, n, :]
            y = self.history[2, n, :]
            z = self.history[3, n, :]

            u = self.history[4, n, :]
            v = self.history[5, n, :]
            w = self.history[6, n, :]

            self.ax1_1.plot(z, x)
            self.ax1_2.plot(z, y)

            self.ax2_1.plot(t, u)
            self.ax2_2.plot(t, v)
            self.ax2_3.plot(t, w)

            self.ax3_1.scatter(x[0], y[0])
            self.ax3_2.scatter(x[-1], y[-1])
        
        # Axis Aesthetics
        self.ax1_1.set_title("XZ Space")
        self.ax1_2.set_title("YZ Space")

        self.ax2_1.set_title("XT Space")
        self.ax2_2.set_title("YT Space")
        self.ax2_3.set_title("ZT Space")

        self.ax3_1.set_title("Initial Position")
        self.ax3_2.set_title("Final Position")

        # Random Actions
        """
        self.ax1_1.axvline(x=l, color='r', linestyle='--')
        self.ax1_2.axvline(x=l, color='r', linestyle='--')

        self.ax1_1.axhline(y=0, color='g', linestyle='--')
        self.ax1_2.axhline(y=0, color='g', linestyle='--')

        self.ax1_1.axhline(y=r, color='g', linestyle='--')
        self.ax1_2.axhline(y=r, color='g', linestyle='--')

        self.ax1_1.axhline(y=-r, color='g', linestyle='--')
        self.ax1_2.axhline(y=-r, color='g', linestyle='--')
        """

        # Create results directory, if missing
        os.makedirs("results", exist_ok=True)
        os.makedirs("results/" + self.date, exist_ok=True)

        if save_fig:
            # Create results directory, if missing
            os.makedirs("results", exist_ok=True)
            os.makedirs("results/" + self.date, exist_ok=True)

            # Save the figures
            self.fig1.savefig("results/"+ self.date + "/SpaceCurve " + groupCode + " "  + self.date + " " + self.hour + ".svg", transparent=True)
            self.fig2.savefig("results/"+ self.date + "/VelocityEvolution " + groupCode + " "  + self.date + " " + self.hour + ".svg", transparent=True)
            self.fig3.savefig("results/"+ self.date + "/HeatMap " + groupCode + " "  + self.date + " " + self.hour + ".svg", transparent=True)

        plt.tight_layout()
        plt.show()
    
    def calibrate(self, gizmo):
        
        # Check for loss[rC] and zone[zC]
        r = gizmo.radialCheck(self.state[1, :], self.state[2, :]) # true if IN r_o
        z = gizmo.zoneCheck(self.state[3, :]) # give you the zone index + 1

        # If the ion is lost, set that in L-vector
        L = r * z 
        print(L, z, r)

        # Ensure zero velocity on loss or detection
        C = np.isin(L, [0, gizmo.NUMBER_OF_MODELS + 1])

        #print(L, r, z, C)

        self.state[0, :] = L
        self.state[4:7, :] *= ~C[None, :]

    def turn(self, gizmo):

        T = len(self.time)

        for t in range(T - 1):

            # Run forward calculations at each time point - the last
            self.history[:, :, t] = self.state
            self.solver.predict(gizmo, self.time[t])

            # Check the zone of each ion in the species set
            self.calibrate(gizmo)
        
        # Add forward state to the final time instant
        self.history[:, :, -1] = self.state

class Simulator:

    def __init__(self):

        # State storage attribute
        self.group = {}
        self.solver = None
        self.date = None
        self.hour = None
        
        # Create gizmo attribute
        self.gizmo = None

        # Create the Matplotlib figures
        self.fig1, (self.ax1_1, self.ax1_2) = plt.subplots(2, 1)
        self.fig2, (self.ax2_1, self.ax2_2, self.ax2_3) = plt.subplots(3, 1)
        self.fig3, (self.ax3_1, self.ax3_2) = plt.subplots(2, 1)
    
    def groupUp(self, name):


        filename_with_ext = name.split("/")[-1]

        self.name = filename_with_ext.split(".")[0]

        with open(name) as f:
            
            i = 0
            rows = csv.reader(f, delimiter=",")

            for row in rows:
                i += 1 
                
                if i==1:
                    pass
                else:
                    tag = row[0]

                    # Initial Definition of Group Member
                    self.group[tag] = State()
                    self.group[tag].tag = tag

                    # Loading ion species parameters since we are here anyway
                    self.group[tag].number = int(row[1])
                    self.group[tag].mass = float(row[2]) * 1.66e-27
                    self.group[tag].charge = float(row[3]) * 1.6e-19
                    self.group[tag].sprX = float(row[4]) * 1e-3
                    self.group[tag].sprY = float(row[5]) * 1e-3
                    self.group[tag].sprVX = float(row[6]) * 1e-3
                    self.group[tag].sprVY = float(row[7]) * 1e-3
                    self.group[tag].sprVZ = float(row[8]) * 1e3
                    self.group[tag].injSpeed = float(row[9]) * 1e3
                    self.group[tag].iniPhase = float(row[10])

                    self.group[tag].solver = self.solver(self.group[tag])
    
    def clock(self):
        
        # Keep track of current time
        pureDate = datetime.now()
        s = pureDate.strftime("%Y-%m-%d %H_%M_%S")
        
        self.date = s[0:10]
        self.hour = s[11:19]

        return self.date, self.hour

    def round(self):

        gizmo = self.gizmo

        self.clock()



        # Loop through every single state
        for code in self.group:

            state = self.group[code]
            state.date, state.hour = self.date, self.hour

            state.h = gizmo.period() * 1e-2
            T = int(1.05 * (gizmo._dst[-1] / state.injSpeed) // state.h)

            gizmo.build(state, state.number)
            state.build(T)

            state.turn(gizmo)
        
        self.l = gizmo._dst[-1]
        self.r = gizmo.inscRadius
    
    def save(self):
        
        # Loop through every single state
        for code in self.group:

            state = self.group[code]

            state.save(groupCode=code)
    
    def display(self, save_fig=True):

        l = self.l
        r = self.r

        # Loop through every single state
        for code in self.group:

            state = self.group[code]
            
            N = state.number
            t = state.time

            for n in range(N):

            # Assignments
                
                x = state.history[1, n, :]
                y = state.history[2, n, :]
                z = state.history[3, n, :]

                u = state.history[4, n, :]
                v = state.history[5, n, :]
                w = state.history[6, n, :]

                self.ax1_1.plot(z, x)
                self.ax1_2.plot(z, y)

                self.ax2_1.plot(t, u)
                self.ax2_2.plot(t, v)
                self.ax2_3.plot(t, w)

                self.ax3_1.scatter(x[0], y[0])
                self.ax3_2.scatter(x[-1], y[-1])
        
        # Axis Aesthetics
        self.ax1_1.set_title("XZ Space")
        self.ax1_2.set_title("YZ Space")

        self.ax2_1.set_title("XT Space")
        self.ax2_2.set_title("YT Space")
        self.ax2_3.set_title("ZT Space")

        self.ax3_1.set_title("Initial Position")
        self.ax3_2.set_title("Final Position")

        # Random Actions
        self.ax1_1.axvline(x=l, color='r', linestyle='--')
        self.ax1_2.axvline(x=l, color='r', linestyle='--')

        self.ax1_1.axhline(y=0, color='g', linestyle='--')
        self.ax1_2.axhline(y=0, color='g', linestyle='--')

        self.ax1_1.axhline(y=r, color='g', linestyle='--')
        self.ax1_2.axhline(y=r, color='g', linestyle='--')

        self.ax1_1.axhline(y=-r, color='g', linestyle='--')
        self.ax1_2.axhline(y=-r, color='g', linestyle='--')

        # Create results directory, if missing
        os.makedirs("results", exist_ok=True)
        os.makedirs("results/" + self.date, exist_ok=True)

        if save_fig:
            # Create results directory, if missing
            os.makedirs("results", exist_ok=True)
            os.makedirs("results/" + self.date, exist_ok=True)

            # Save the figures
            self.fig1.savefig("results/"+ self.date + "/SpaceCurve "+ self.name + " "  + self.date + " " + self.hour + ".svg", transparent=True)
            self.fig2.savefig("results/"+ self.date + "/VelocityEvolution " + self.name + " "  + self.date + " " + self.hour + ".svg", transparent=True)
            self.fig3.savefig("results/"+ self.date + "/HeatMap " + self.name + " "  + self.date + " " + self.hour + ".svg", transparent=True)

        plt.tight_layout()
        plt.show()



