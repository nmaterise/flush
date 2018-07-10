#%%
import matplotlib.pyplot as pt
import numpy as np

filename = "C:/Users/Bebotron/Documents/Research/flush/three_qubit/"

[t, F] = np.loadtxt(filename + "output_20_40.dat")

pt.plot(t, F)
pt.show()