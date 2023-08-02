# %%
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('NodeDisp.out')
plt.plot(data[:, 0], data[:, 1])
plt.plot(data[:, 0], data[:, 3])
plt.plot(data[:, 0], data[:, 5])
plt.xlabel('Time (s)')
plt.ylabel('Z Displacement (m)')
plt.legend(labels=['1', '2', '3'])
plt.ylim(-0.0006, 0.0006)
plt.grid(linestyle='--')
plt.title('PML single core ')
plt.show()

# %%
# data = np.loadtxt('NodeDispx.out')
# plt.plot(data[:, 0], data[:, 1])
# plt.plot(data[:, 0], data[:, 3])
# plt.plot(data[:, 0], data[:, 5])
# plt.xlabel('Time (s)')
# plt.ylabel('X Displacement (m)')
# plt.legend(labels=['1', '2', '3'])
# plt.grid(linestyle='--')
# plt.ylim(-0.0006, 0.0006)
# plt.title('With PML single core ')
# %%
