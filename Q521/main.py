import numpy as np
from matplotlib import pyplot as plt
import math
from qmf import QMF
from history import History
from corral import Corral
from law import Matthieu_Equation
from solver import Euler


# Creating the QMF and History Objects
qmf = QMF("filter1.csv", "McIntosh Filter")
qmf.load_csv()


f = qmf.pureFrequency
v = qmf.injSpeed
l = 5e-3 #qmf.span
T = 1/(100 * f)
print( 100 * (l/v)//T)
N_STEPS = int(100 * (l/v)//T) + 50
N_NODES = N_STEPS + 1


ionA = History("ionspec1.csv", "Sodium")
ionA.load_csv()
ionA.setup(N_NODES, qmf)

poleField = Corral(Euler, Matthieu_Equation, T, N_STEPS, qmf, ionA)
poleField.lockOn()
t, x, v, a, q = poleField.main()
print(t[0, :])





"""
Data Visualization
"""


#"""

t, x, v, a, q = poleField.main()

hstTime = t[0, :]
x1 = x[0, :, 0]
x2 = x[1, :, 0]
r = np.sqrt((x1 ** 2) + (x2 ** 2))


# Create Figures

fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
ax1.plot(hstTime, x1)
ax1.set_title("XZ Space")
ax2.plot(hstTime, x2)
ax2.set_title("YZ Space")
ax3.plot(hstTime, r)
ax3.set_title("RT Space")

plt.tight_layout()
plt.show()
#"""
