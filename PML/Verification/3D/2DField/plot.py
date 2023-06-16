# %%
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('NodeDisp.out')
plt.plot(data[:, 0], data[:, 2])
plt.plot(data[:, 0], data[:, 6])
# plt.plot(data[:, 0], data[:, 5])

# %%
