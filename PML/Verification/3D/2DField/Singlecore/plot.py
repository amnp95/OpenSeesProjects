# %%
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('NodeDisp.out')
plt.plot(data[:, 0], data[:, 1], label = 1)
# plt.plot(data[:, 0], data[:, 3], label = 2)
plt.plot(data[:, 0], data[:, 5], label = 3)
plt.legend()
# %%
