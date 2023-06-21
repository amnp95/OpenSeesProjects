# %%
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('NodeDisp.out')
plt.plot(data[:, 0], data[:, 3])
plt.plot(data[:, 0], data[:, 1])
# plt.plot(data[:, 0], data[:, 5])

# %%
