# %%
import pyvista as pv
import numpy as np
import os 

# change the working directory to the file location
if os.path.dirname(__file__) != '':
    os.chdir(os.path.dirname(__file__))

case = "PML"

pileforce = np.loadtxt(f"./results/{case}/BeamForce0.out") 


# %%
fx = np.concatenate((pileforce[:,0::12],-pileforce[:,-6].reshape(-1,1)), axis= 1)
fy = np.concatenate((pileforce[:,1::12],-pileforce[:,-5].reshape(-1,1)), axis= 1)
fz = np.concatenate((pileforce[:,2::12],-pileforce[:,-4].reshape(-1,1)), axis= 1)
Mx = np.concatenate((pileforce[:,3::12],-pileforce[:,-3].reshape(-1,1)), axis= 1)
My = np.concatenate((pileforce[:,4::12],-pileforce[:,-2].reshape(-1,1)), axis= 1)
Mz = np.concatenate((pileforce[:,5::12],-pileforce[:,-1].reshape(-1,1)), axis= 1)

# %%

time = np.loadtxt(f"./results/{case}/Time.out")[:,0]

import matplotlib.pyplot as plt
depth = np.loadtxt(f"/home/amnp95/Projects/OpenSeesProjects/PML_DRM_Pile/BeamEqStruct/results/beamNodes.dat")[:,-1]
plt.plot(My.T,depth, color = 'blue')
plt.show()

# change M from N.m to kN.m
Mx = Mx/1e3
My = My/1e3
Mz = Mz/1e3

# change F from N to kN
fx = fx/1e3
fy = fy/1e3
fz = fz/1e3


# %%
# plot the envelope of the moment
plt.plot(np.max(My, axis = 0),depth, color = 'blue',alpha = 0.2)
plt.plot(np.min(My, axis = 0),depth, color = 'blue',alpha = 0.2)
# fill the envelope
plt.fill_betweenx(depth, np.max(My, axis = 0), np.min(My, axis = 0), color = 'blue', alpha = 0.2)

plt.xlabel('Moment (kN.m)')
plt.ylabel('Depth (m)')

# %%
# create an animation of the moment
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots()
xdata, ydata = [], []
ln, = ax.plot([], [], '-',color = 'blue',alpha = 0.2)
ln2, = ax.plot([], [], 'o-',color = 'red',alpha = 1.0, linewidth = 1, markersize = 3.5)

def init():
    ax.set_xlim(-5000,5000)
    ax.set_ylim(0, 3.5)
    ax.set_xlabel('Moment (kN.m)')
    ax.set_ylabel('Depth (m)')
    ax.grid(linestyle='--')
    return ln,

def update(frame):
    xdata.append(My[frame,:])
    ydata.append(depth)
    ln.set_data(xdata, ydata)
    # clean ln2
    ln2.set_data([],[])
    ln2.set_data(My[frame,:],depth)
    return ln,

ani = FuncAnimation(fig, update, frames=np.arange(300, 1000), init_func=init, blit=True)
ani.save('test2.mp4',  
          writer = 'ffmpeg', fps = 30) 

# %%
# ==================================================================
# beam mesh
# ==================================================================
beamnodes = np.loadtxt("results/beamNodes.dat")
beamcells = np.loadtxt("results/beamElems.dat", dtype=int)
beamnodes = beamnodes[:,1:]
beamcells[:,0] = 2
beamcells[:,1:] = beamcells[:,1:] - 1
celltypes  = np.ones(beamcells.shape[0],dtype= int) * pv.CellType.LINE
beam = pv.UnstructuredGrid(beamcells.tolist(),celltypes.tolist(),beamnodes.tolist())


# %%
# ==================================================================
# create interface mesh
# ==================================================================
data = np.loadtxt("results/interfaceInfo.dat")
interface = pv.PolyData(data[:,1:4])


