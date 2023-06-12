
# %%
import numpy as np
import matplotlib.pyplot as plt



data = np.loadtxt("node_disp.out")
plt.plot(data[:,0], data[:,2])
data = np.loadtxt("node_disp2.out")
plt.plot(data[:,0], data[:,1])
# %%
data