# %%
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('NodeDisp.out')
plt.plot(data[:, 0], data[:, 1], label = "node 1")
# plt.plot(data[:, 0], data[:, 2], label = "node 2")
# plt.plot(data[:, 0], data[:, 3], label = "Footing")

# analytical = np.loadtxt('analytical.csv', delimiter = ',')
# plt.plot(analytical[:, 0], analytical[:, 1]*1e-3, label = "Analytical")

plt.xlabel('Time (s)')
plt.ylabel(' Vertical Displacement (m)')
plt.legend()
# plt.title("Footing Displacment using PML")
# plt.title("Footing Displacment Using Fixed Boundaries at base and side PMLs")
# plt.title("Footing Displacment Using Fixed Boundaries (f=8Hz))")
# plt.title("Footing Displacment Using PML Boundaries")


plt.grid(linestyle = '--')
# plt.ylim(-0.0045, 0.007)
plt.show()








# %%
import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt('NodeDisp.out')
plt.plot(data[:, 0], data[:, 3], label = "Footing")
data = np.loadtxt('Luco_vv_025_tipo_de_carga_11_duracion_0.4RichterHF.out',delimiter=',')
plt.plot(data[:, 0], data[:, 1],label = "analytical")
plt.xlabel('Time (s)')
plt.ylabel(' Vertical Displacement (m)')
plt.title("Footing Displacment using PML")
plt.xlim(0, 1.0)
plt.grid(linestyle = '--')
plt.legend()
plt.show()
# %%
