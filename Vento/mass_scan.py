import time
from sim.state import State, Simulator
from sim.gizmo import Gizmo, Ideal_Zone, Detector_Zone
from sim.solver import Euler, RK4


sim = Simulator()
sim.solver = RK4
sim.groupUp(name = "sim/ion/Set A.csv")

machine = Gizmo()
machine.tag = "Reis Filter"
machine.load("sim/filter/Standard Series.csv")

A = Ideal_Zone(machine, 20)
B = Ideal_Zone(machine, 600)
C = Ideal_Zone(machine, 20)
D = Detector_Zone(machine)

machine._val = [A, B, C, D]
machine._dst = [0.02, 0.62, 0.64]

sim.gizmo = machine
sim.round()

sim.save()
sim.display()

print(sim.group)

print(sim.group['Titanium-47.5'].mass)