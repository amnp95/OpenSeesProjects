# %%
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('NodeDisp.out')
plt.plot(data[:, 0], data[:, 1], label = "node 1")
plt.plot(data[:, 0], data[:, 3], label = "node 2")
plt.plot(data[:, 0], data[:, 5], label = "node 3")
plt.xlabel('Time (s)')
plt.ylabel(' Vertical Displacement (m)')
plt.legend()
plt.title("2D model with PML")
plt.grid(linestyle = '--')
plt.show()
# %%

