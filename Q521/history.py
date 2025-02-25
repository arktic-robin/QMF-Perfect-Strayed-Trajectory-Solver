import math
import numpy as np
import csv
from copy import deepcopy
import random as rng
from qmf import QMF

class History:
    """
    This class is intended to store the parameters and coordinates
    of a specific grouping of ions through time or phase
    """
    
    def __init__(self, name, tag, param=[]):
        # Input Parameters
        self.name = name
        self.param = None
        self.tag = tag
        if name == None:
            self.tag = param[0]
            self.number = param[1]
            self.mass = param[2] * 1.66e-27
            self.charge = param[3] * 1.6e-19
            self.spX = param[4] * 1e-3
            self.spY = param[5] * 1e-3
            self.spVX = param[6] * 1e-3
            self.spVY = param[7] * 1e-3

    def setup(self, nodes, filter):
        # Array Initialization
        self.phaseTime = np.zeros([2, nodes])
        self.pos = np.zeros([3, nodes, self.number])
        self.vel = np.zeros([3, nodes, self.number])
        self.acl = np.zeros([3, nodes, self.number])
        self.gen = np.zeros([2, nodes, self.number])

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

            self.acl[0, 0, i] = 0
            self.acl[1, 0, i] = 0
            self.acl[2, 0, i] = 0

            self.gen[0, 0, i] = math.sqrt( (del_x ** 2) + (del_y ** 2)) / filter.inscRadius
            self.gen[1, 0, i] = math.sqrt( (del_vx ** 2) + (del_vy ** 2)) / filter.refVelocity

    
    def load_csv(self):
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
        clone = self.replicate()
        clone.phaseTime = clone.phaseTime[:, j]
        clone.pos = clone.pos[:, j, :]
        clone.vel = clone.vel[:, j, :]
        clone.acl = clone.acl[:, j, :]
        clone.gen = clone.gen[:, j, :]
        return clone
    
    def replicate(self):
        clone = deepcopy(self)
        return clone
    
    def step(self, j, dt, w):
        self.phaseTime[0, j] = self.phaseTime[0, j - 1] + dt
        self.phaseTime[1, j] = self.phaseTime[1, j - 1] + w * dt
