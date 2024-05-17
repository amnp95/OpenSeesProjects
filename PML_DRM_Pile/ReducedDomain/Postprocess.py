# %%
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt 
import os
import tqdm
# change the path to the directory where the file is located
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# ==================================================================
# load nodes and elements dataframes
# ==================================================================
mesh = pv.read("results/mesh.vtk")
info = {
    "RegularDomain": 1,
    "DRMDomain": 2,
    "PMLDomain": 3,
}
# filter cells with Domin tag 1
indices = mesh['Domain'] == 1
grid = mesh.extract_cells(indices)
cleangrid = grid.copy()
cleangrid.clear_data()
# %%
# ==================================================================
# load displacement
# ==================================================================
# typename = "PML"
typename = "FIXED"
Dir = f"results/{typename}"
griddisp   = np.loadtxt(f"{Dir}/NodeDisp0.out")
gridaccel  = np.loadtxt(f"{Dir}/NodeAccl0.out")
time       = np.loadtxt(f"{Dir}/Time.out")[:,0]
# %%
dt        = time[1] - time[0]
frameRate = 1/dt
griddispx = griddisp[:,0::3]
griddispy = griddisp[:,1::3]
griddispz = griddisp[:,2::3]
dispfactor = 50
# %%
#  ==================================================================
#  create  directory to save the images
#  ==================================================================

direcotry = f"./movie/Pics{typename}"
if not os.path.exists(direcotry):
    if not os.path.exists("./movie"):
        os.mkdir("./movie")
    os.mkdir(direcotry)

# ==================================================================
# slice the mesh
# ==================================================================

slicer = pv.Cube(center=[0,-0.25,0], x_length=150, y_length=0.1, z_length=150)
slicedmesh = cleangrid.copy()
slicedmesh = slicedmesh.clip_box(slicer,invert=False,crinkle=True)
slicepointindexes = slicedmesh["vtkOriginalPointIds"]
slicecelliindexes = slicedmesh["vtkOriginalCellIds"]

originalcoords     = slicedmesh.points
name  =  "acceleration"
# array =  np.abs(gridaccel[:,0::3])
array =  gridaccel[:,0::3]
mx    =  array.max()
mn    =  array.min()
mn,mx =  np.sort(np.abs([mx,mn]))
array =  array[:,slicepointindexes]
# frameRate = 50
mx = 10
mn = -5

# %%
# ==================================================================
# create the movie
# ==================================================================
# empty the directory
files = os.listdir(direcotry)
for file in files:
    os.remove(f"{direcotry}/{file}")

pl = pv.Plotter(off_screen=True)     


# create progress bar
step  = 1
start = 330
# end   = griddisp.shape[0]
end   = 1000
total = end - start
total = total//step
bar = tqdm.tqdm(total=total)


for i in range(start,end,step):
    bar.update(1)
    slicedmesh.points = originalcoords     + dispfactor * griddisp[i,:].reshape(-1,3)[slicepointindexes]
    slicedmesh.point_data[f"{name}"] = array[i,:]
    pl.add_mesh(slicedmesh,scalars=f"{name}", show_edges=True, cmap="rainbow", clim=[mn,mx])
    pl.camera_position = 'xz'
    pl.save_graphic(f"{direcotry}/{i:04d}.svg")
    pl.clear()
pl.close()
bar.close()

# %%
direcotry = f"./movie/Pics{typename}"
frameRate = 30
os.system(f"ffmpeg -y -framerate {frameRate} -pattern_type glob -i '{direcotry}/*.svg' -c:v libx264 -pix_fmt yuv420p tmp.mp4")

# %%
