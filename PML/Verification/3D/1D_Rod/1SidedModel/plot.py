
# %%
import numpy as np
import matplotlib.pyplot as plt

# %%
data = np.loadtxt("NodeDispHead.out")
plt.plot(data[:,0], data[:,2])
data = np.loadtxt("NodeDispEnd.out")
plt.plot(data[:,0], data[:,9])
plt.grid(linestyle='--')
# %%

# %%

data = np.loadtxt("node_disp1.out")
plt.plot(data[:,0], data[:,2])
data = np.loadtxt("node_disp2.out")
plt.plot(data[:,0], data[:,1])
# %%
