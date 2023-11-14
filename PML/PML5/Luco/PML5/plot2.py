# %%
import numpy as np
import matplotlib.pyplot as plt



Hf = np.loadtxt("Luco_vv_025_tipo_de_carga_11_duracion_0.4RichterHF.out", delimiter=",")
Lf = np.loadtxt("Luco_vv_025_tipo_de_carga_11_duracion_1RichterLF.out", delimiter=",")
# %%
plt.plot(Hf[:,0],Hf[:,1],label="High frequncy loading")
# plt.plot(Lf[:,0],Lf[:,1],label="Low frequncy loading")
plt.xlim(right= 1.0)
plt.xlim(left= 0.0)
plt.title ("High frequncy loading")
plt.xlabel("Time [s]")
plt.ylabel("Force  [kN]")
plt.grid(linestyle='--')
plt.ylim(bottom= -0.005)
plt.ylim(top= 0.007)


# %%
plt.plot(Lf[:,0],Lf[:,1],label="Low frequncy loading")
plt.xlim(right= 1.0)
plt.xlim(left= 0.0)
plt.title ("Low frequncy loading")
plt.xlabel("Time [s]")
plt.ylabel("Force  [kN]")
plt.grid(linestyle='--')
plt.ylim(bottom= -0.005)
plt.ylim(top= 0.007)
plt.show()
# %%

data = np.loadtxt("NodeDisp0.out")
plt.plot(data[:,0],data[:,1],label="PML 5")
data = np.loadtxt("/home/amin/Projects/Github/OpenSeesProjects/PML/PML5/Luco/PML2D/NodeDisp0.out")
plt.plot(data[:,0],data[:,1],label="PML 2D")

Hf = np.loadtxt("Luco_vv_025_tipo_de_carga_11_duracion_0.4RichterHF.out", delimiter=",")
Lf = np.loadtxt("Luco_vv_025_tipo_de_carga_11_duracion_1RichterLF.out", delimiter=",")
plt.plot(Hf[:,0],Hf[:,1],label="Analytical")
plt.grid(linestyle='--')
plt.title ("Displacement for High frequncy loading at node 1") 
plt.xlabel("Time [s]")
plt.ylabel("Displacement  [m]")
plt.legend()
plt.xlim(right= 1.0)
plt.xlim(left= 0.0)
# %%
data = np.loadtxt("/home/amin/Projects/Github/OpenSeesProjects/PML/PML5/Luco/PML5/RichterHF.dat")
plt.plot(np.arange(0, data.shape[0]*0.001, 0.001),data)
plt.xlim(right= 1.0)
plt.xlim(left= 0.0)
plt.title("High Frequncy Loading")
plt.xlabel("Time [s]")
plt.ylabel("Force  [kN]")
plt.grid(linestyle='--')

# %%
