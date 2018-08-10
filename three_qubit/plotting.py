#%%
import matplotlib.pyplot as pt
import numpy as np

filename = "C:/Users/Bebotron/Documents/Research/flush/three_qubit/"

[t, F1, F2, F3] = np.loadtxt(filename + "sbatchFiles/outFiles/yes_coupling.dat")

pt.plot(t, F1)
# pt.plot(t, F2)
# pt.plot(t, F3)
pt.show()