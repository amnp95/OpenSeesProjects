# %%
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt 
import os
# change the path to the directory where the file is located
os.chdir(os.path.dirname(os.path.abspath(__file__)))


# ==================================================================
# functions
# ==================================================================
# helper function to load data
def loadData3D(baseName,directory,nodeTagsfile,eleTagsfile,datashape):
    files = os.listdir(directory)
    basefiles = sorted([file for file in files if baseName in file])
    nodetags  = sorted([file for file in files if nodeTagsfile in file])
    eletags   = sorted([file for file in files if eleTagsfile in file])

    nodetags = [np.loadtxt(f"{directory}/{file}",dtype=int)-1 for file in nodetags]
    eletags  = [np.loadtxt(f"{directory}/{file}",dtype=int)-1 for file in eletags]
    data     = np.empty(datashape,dtype=np.float64)

    for i in range(len(basefiles)):
        tmpdata = np.loadtxt(f"{directory}/{basefiles[i]}")[:,1:]
        indices = np.array([3 * i + np.array([0, 1, 2]) for i in nodetags[i]]).flatten()
        indices = indices.tolist()

        # check that the indicies length is the same as the data length
        if len(indices) != tmpdata.shape[1]:
            raise ValueError(f"the indices length {len(indices)} is different from the data length {tmpdata.shape[1]}")
        
        # assign the data to the correct location
        data[:,indices] = tmpdata
    return data 

# ==================================================================
# load nodes and elements dataframes
# ==================================================================
mesh = pv.read("results/mesh.vtk")
info = {"RegularDomain": 1, "DRMDomain": 2, "PMLDomain": 3}

# filter cells with Domin tag 1
indices = mesh['Domain'] == 1
grid = mesh.extract_cells(indices)
cleangrid = grid.copy()
cleangrid.clear_data()
# %%
tmp  = mesh.extract_cells(np.where(mesh.cell_data['partitioned']==2)[0])
# tmp.plot(show_edges=True,color="blue")
tmp.n_points
direcotry = "results"
tmpdata0 = np.loadtxt(f"{direcotry}/NodeDisp2.out")
tmpdata1 = np.loadtxt(f"{direcotry}/NodeDispPML2.out")
index = tmp.find_closest_point([0,0,0])
plt.plot(tmpdata1[:,0],tmpdata0[:,index*3+1],"b",label="Fixed boundary")
plt.plot(tmpdata0[:,0],tmpdata1[:,index*3+1],"r",label="PML boundary")
plt.legend(loc="upper right",framealpha=1.0)
plt.xlabel("Time [s]")
plt.ylabel("Displacement [m]")
plt.xlim([0,4])
plt.ylim([-0.015,0.015])

# %%
# ==================================================================
# load data
# ==================================================================
time    = np.loadtxt(f"results/NodeAcclPML0.out")[:,0]
npoints = grid.n_points
steps   = time.shape[0]
indices = np.array([3 * i + np.array([0, 1, 2]) for i in grid["vtkOriginalPointIds"]]).flatten().tolist()

data    = loadData3D("NodeAcclPML","results","nodeTags","elementTags",datashape=(steps,npoints*3))
disp    = loadData3D("NodeDispPML","results","nodeTags","elementTags",datashape=(steps,npoints*3))
data    = data[:,indices]
disp    = disp[:,indices]
data    = data[:,0::3]

# %%
dt        = time[1] - time[0]
frameRate = 1/dt
griddispx = disp[:,0::3]
griddispy = disp[:,1::3]
griddispz = disp[:,2::3]
dispfactor = 50
# %%
slicer = pv.Cube(center=[0,-0.25,0], x_length=150, y_length=0.24, z_length=150)
slicedmesh = cleangrid.copy()
slicedmesh = slicedmesh.clip_box(slicer,invert=False,crinkle=True)
slicepointindexes = slicedmesh["vtkOriginalPointIds"]
slicecelliindexes = slicedmesh["vtkOriginalCellIds"]
# %%

direcotry = "./movie/Pics"
if not os.path.exists(direcotry):
    if not os.path.exists("./movie"):
        os.mkdir("./movie")
    os.mkdir(direcotry)
originalcoords     = slicedmesh.points
data  =  np.abs(data)
data  = data[:,slicepointindexes]
# %%

files = os.listdir(direcotry)
for file in files:
    os.remove(f"{direcotry}/{file}")
name = "acceleration"
pl = pv.Plotter(off_screen=True)     
mx = 1.0
mn = 0.0
frameRate = 60
# for i in range(500,800,1):
for i in range(0,steps,1):
    slicedmesh.points = originalcoords     + dispfactor * disp[i,:].reshape(-1,3)[slicepointindexes]
    slicedmesh.point_data[f"{name}"] = data[i,:]
    pl.add_mesh(slicedmesh,scalars=f"{name}", show_edges=True, cmap="coolwarm", clim=[mn,mx])
    pl.camera_position = 'xz'
    pl.save_graphic(f"{direcotry}/{i:04d}.svg")
    pl.clear()
pl.close()

# %%
import os 
direcotry = "./movie/Pics"
frameRate = 60
os.system(f"ffmpeg -y -framerate {frameRate} -pattern_type glob -i '{direcotry}/*.svg' -c:v libx264 -pix_fmt yuv420p tmp.mp4")



# %% 
tmp = mesh["partitioned"] == 2
grid = mesh.extract_cells(tmp)

i = grid.find_closest_point([0,0,3])

data = np.loadtxt("results/NodeDispPML2.out")
time = data[:,0]
data = data[:,1::3]

plt.plot(time,data[:,i])

data = np.loadtxt("results/NodeDisp2.out")
time = data[:,0]
data = data[:,1::3]
plt.plot(time,data[:,i])















# %%
