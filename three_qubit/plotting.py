#%%
import matplotlib.pyplot as pt
import numpy as np

filename = "C:/Users/Bebotron/Documents/Research/flush/three_qubit/"

[t, F] = np.loadtxt(filename + "outputF2_20.dat")

pt.plot(t, F)
pt.show()