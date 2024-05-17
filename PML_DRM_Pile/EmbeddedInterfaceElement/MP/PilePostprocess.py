# %%
import pyvista as pv
import numpy as np
import os 

# change the working directory to the file location
if os.path.dirname(__file__) != '':
    os.chdir(os.path.dirname(__file__))

case = "PML"

pileforce = np.loadtxt(f"./results/BeamForce{case}0.out") 


# %%
fx = np.concatenate((pileforce[:,1::12],-pileforce[:,-6].reshape(-1,1)), axis= 1)
fy = np.concatenate((pileforce[:,2::12],-pileforce[:,-5].reshape(-1,1)), axis= 1)
fz = np.concatenate((pileforce[:,3::12],-pileforce[:,-4].reshape(-1,1)), axis= 1)
Mx = np.concatenate((pileforce[:,4::12],-pileforce[:,-3].reshape(-1,1)), axis= 1)
My = np.concatenate((pileforce[:,5::12],-pileforce[:,-2].reshape(-1,1)), axis= 1)
Mz = np.concatenate((pileforce[:,6::12],-pileforce[:,-1].reshape(-1,1)), axis= 1)

# %%

time = pileforce[:,0]

import matplotlib.pyplot as plt
depth = np.loadtxt(f"./results/beamNodes.dat")[:,-1]
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
ln, = ax.plot([], [], '-',color = 'blue',alpha = 0.2, linewidth = 2)
ln2, = ax.plot([], [], 'o-',color = 'red',alpha = 1.0, linewidth = 1, markersize = 3.5)
xmax = np.max(My)*1.1
xmin = np.min(My)*1.1
ymax = np.max(depth)+0.5
ymin = np.min(depth)-0.5
ax.hlines(0,xmin*2,xmax*2,linestyle='-',color = 'black',linewidth = 2)
ax.vlines(0,ymin*2,ymax*2,linestyle='--',color = 'black',linewidth = 0.5)  
def init():
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel('Moment (kN.m)')
    ax.set_ylabel('Depth (m)')
    ax.grid(linestyle='--')
    return ln,

def update(frame):
    xdata.append(My[frame,:])
    ydata.append(depth)
    ln.set_data(xdata, ydata)
    # clean ln2
    ln2.set_data(My[frame,:],depth)
    return ln,

ani = FuncAnimation(fig, update, frames=np.arange(300, 1000), init_func=init, blit=True)
ani.save(f'Pilemoment{case}.mp4',  
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
beamdisp = np.loadtxt(f"./results/BeamDisp{case}0.out")[:,1:]
beamoriginalcoord = beam.points
# %%
# ==================================================================
# create interface mesh
# ==================================================================
data = np.loadtxt("results/interfaceInfo.dat")
x = data[:,1].reshape(-1,8)
y = data[:,2].reshape(-1,8)
z = data[:,3].reshape(-1,8)

# add the first 8 points in each row to the end of the row
x = np.concatenate((x,x[:,0].reshape(-1,1)), axis = 1)
y = np.concatenate((y,y[:,0].reshape(-1,1)), axis = 1)
z = np.concatenate((z,z[:,0].reshape(-1,1)), axis = 1)
interface = pv.StructuredGrid(x,y,z).merge()

# Create a dictionary to map each row of the original matrix to its index
row_to_index = {tuple(row): idx for idx, row in enumerate(data[:,1:4])}
interface_mapping = [row_to_index[tuple(row)] for row in interface.points]
# for each element add the (3*i, 3*i+1, 3*i+2) to the mapping consecutively
interface_mapping = [3*i+j for i in interface_mapping for j in range(3)]

interfacedisp = np.loadtxt(f"./results/Interfacepoints{case}0.out")[:,1:]
interfacedisp = interfacedisp[:,interface_mapping]

interfaceoriginalcoord = interface.points


# ==================================================================
# create alumped mass
# ==================================================================
mass = beam.points[-1,:]
mass = pv.PolyData(mass.reshape(1,-1))
massdisp = beamdisp[:,[-3,-2,-1]]
massoriginalcoord = mass.points
#%%
pl = pv.Plotter()
pl.add_mesh(beam, color = 'black', line_width = 5)
pl.add_mesh(interface, show_edges=True, line_width = 1, show_vertices=True, vertex_color='red',use_transparency=True,opacity=0.5, render_points_as_spheres = True)
pl.add_mesh(mass, color = 'blue', point_size = 50, render_points_as_spheres = True)
pl.camera_position = [(10.382977869678236, 10.881033743089947, 6.407282603898192),
                      (0.0, 0.0, -1.5),
                      (-0.30520760684173853, -0.35182912615338785, 0.8849093641249831)]
pl.show()
# %%
dispfactor = 50
pl = pv.Plotter(notebook=False, off_screen=True)
pl.add_mesh(beam, color = 'black', line_width = 5)
pl.add_mesh(interface, show_edges=True, line_width = 1, show_vertices=True, vertex_color='red',use_transparency=True,opacity=0.5, render_points_as_spheres = True)
pl.add_mesh(mass, color = 'blue', point_size = 50, render_points_as_spheres = True)
pl.camera_position = [(10.382977869678236, 10.881033743089947, 6.407282603898192),
                      (0.0, 0.0, -1.5),
                      (-0.30520760684173853, -0.35182912615338785, 0.8849093641249831)]

# pl.open_gif(f"pile{case}.gif",fps=30)
pl.open_movie(f"pile{case}.mp4",framerate=30, quality=10)
for i in range(300,2000):
    beam.points = beamoriginalcoord + beamdisp[i,:].reshape(-1,3)*dispfactor
    interface.points = interfaceoriginalcoord + interfacedisp[i,:].reshape(-1,3)*dispfactor
    mass.points = massoriginalcoord + massdisp[i,:].reshape(-1,3)*dispfactor
    pl.write_frame()

pl.close()

# %%
