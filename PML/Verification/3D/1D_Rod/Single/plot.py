
# %%
import numpy as np
import matplotlib.pyplot as plt

# %%
data = np.loadtxt("NodeDispCentre.out")
plt.plot(data[:,0], data[:,1])
data = np.loadtxt("NodeDispNegative.out")
plt.plot(data[:,0], data[:,1])
data = np.loadtxt("NodeDispPositive.out")
plt.plot(data[:,0], data[:,1])
plt.grid(linestyle='--')
plt.xlabel("Time [s]")
plt.ylabel("x Displacement [m]")
plt.legend(labels=['Center', 'Negative Head', 'Positive Head'])
plt.title("Two Cores With PML in Positive Head and loading on Negative Head")
plt.show()
# %%
