# %%
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('NodeDisp.out')
plt.plot(data[:, 0], data[:, 1],label="without PML")
data = np.loadtxt('NodeDisp_PML.out')
plt.plot(data[:, 0], data[:, 3],label="with PML")


plt.xlabel('Time (s)')
plt.ylabel('Z Displacement (m)')
# plt.ylim(-0.0021, 0.001)
plt.xlim(right=5.1)
plt.grid(linestyle='--')
plt.legend()
plt.savefig('figs/centerNodedisplacment.png')
plt.show()
# %%
