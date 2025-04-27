import time
from sim.state import State
from sim.gizmo import Gizmo, Ideal_Zone, Detector_Zone, Dawson_Entry_Zone, HM_Entry_Zone, Dawson_Exit_Zone
from sim.solver import Euler, RK4

T = 20000
tag = "Titanium-47-A"
name = "sim/ion/Set A.csv"

ion = State()
ion.tag = tag
ion.load(name)
ion.solver = RK4(ion)



machine = Gizmo()
machine.tag = "Reis Filter"
machine.load("sim/filter/Standard Series.csv")

A = HM_Entry_Zone(machine, 20)
B = Ideal_Zone(machine, 300)
C = Dawson_Exit_Zone(machine, 50)
D = Detector_Zone(machine)

machine._val = [A, B, C, D]
l = 0
for val in machine._val[0:-1]:
    l += val.span
    machine._dst.append(l)

print(machine._dst, machine.inscRadius)

S = time.time()

ion.h = machine.period() * 1e-2
T = int(1.05 * (machine._dst[-1] / ion.injSpeed) // ion.h)
print(T)
machine.build(ion, ion.number)
ion.build(T)
ion.clock()
ion.turn(machine)
ion.length = machine._dst[-1]
ion.rad = machine.inscRadius

E = time.time()

print(E - S)

print(machine._dst, machine.inscRadius)
#ion.save()
ion.display()
