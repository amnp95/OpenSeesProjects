
# %%
import numpy as np
import matplotlib.pyplot as plt

# %%
data = np.loadtxt("NodeDispCent.out")
plt.plot(data[:,0], data[:,1])
data = np.loadtxt("NodeDispHead.out")
plt.plot(data[:,0], data[:,1])
data = np.loadtxt("NodeDispEnd.out")
plt.plot(data[:,0], data[:,1])
plt.grid(linestyle='--')
# %%

