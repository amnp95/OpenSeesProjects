# %%
import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt('displacement_select.out')
time = data[:, 0]
ydata = data[:, 2::2]
plt.plot(time, ydata[:, 0], 'r-', label='x')
# %%
