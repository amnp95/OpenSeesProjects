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

# ==================================================================
# create interface mesh
# ==================================================================
data = np.loadtxt("results/interfaceInfo.dat")
interface = pv.PolyData(data[:,1:4])

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

pl = pv.Plotter(shape=(1,2))
pl.subplot(0, 0) 
cube = pv.Cube(center=[0,0,1], x_length=10, y_length=10, z_length=1)
# find cells in the grid inside the cudebox
found = mesh.find_cells_within_bounds(cube.bounds)
# get the cell indexes
mesh["Domain"][found] = 4
found = grid.find_cells_within_bounds(cube.bounds)
grid["Domain"][found] = 4
gridwithoutfound = grid.extract_cells(grid["Domain"] != 4)
found = grid.extract_cells(grid["Domain"] == 4)



pl.add_mesh(gridwithoutfound, color="royalblue", show_edges=True,opacity=1.0,show_scalar_bar=False,cmap=["royalblue","olivedrab","peru","slategrey"],label="Soil")
pl.add_mesh(found, color="grey", line_width=5, opacity=0.5,label="Foundation")
pl.add_mesh(interface, color="green", line_width=5, label="Interface")
pl.add_mesh(beam, color="blue", line_width=5, render_points_as_spheres=True, point_size=8, label="Beam")
pl.add_mesh(beam, style="points", color="red", point_size=7,render_points_as_spheres=True,label="Mass")
pl.add_legend(labels=None, bcolor="white", border=True, size=(0.2, 0.2), name=None, loc='lower right', face='rectangle')
pl.add_axes()

pl.subplot(0, 1)
cube = pv.Cube(center=[0,0.25,0], x_length=100, y_length=0.1, z_length=100)
points = pv.PolyData([[0,0.,3],[0,0,0]])
pl.add_point_labels(points,["       P1 (x=0,y=0,z=3.0)","        P2 (x=0, y=0, z=0)"],font_size=20,point_color="purple",text_color="purple",always_visible=True,shape_color="white",shape_opacity=0.0, render_points_as_spheres=True, point_size=10)
pl.add_mesh(gridwithoutfound.clip_box(cube,invert=False,crinkle=True), color="royalblue", show_edges=True,opacity=1.0,show_scalar_bar=False,cmap=["royalblue","olivedrab","peru","slategrey"],label="Soil")
pl.add_mesh(found.clip_box(cube,invert=False,crinkle=True), color="grey", opacity=1.0,label="Foundation",show_edges=True)
pl.add_mesh(interface, color="green", line_width=5, label="Interface",)
pl.add_mesh(beam, color="blue", line_width=5, render_points_as_spheres=True, point_size=8, label="Beam")
pl.camera_position = 'xz'
pl.add_axes()
pl.add_mesh(points, render_points_as_spheres=True, point_size=10, color="purple",label=" Control points")
pl.add_legend(labels=None, bcolor="white", border=True, size=(0.22, 0.22), name=None, loc='lower right', face='rectangle')
pl.show()



# %%
# ==================================================================
# load displacement
# ==================================================================
typename = "PML"
# typename = "FIXED"
Dir = f"results/{typename}"
griddisp   = np.loadtxt(f"{Dir}/NodeDisp0.out")
beamdisp   = np.loadtxt(f"{Dir}/BeamDisp0.out")
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
beamoriginalcoords = beam.points
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
    beam.points       = beamoriginalcoords + dispfactor * beamdisp[i,:].reshape(-1,3)
    slicedmesh.point_data[f"{name}"] = array[i,:]
    pl.add_mesh(slicedmesh,scalars=f"{name}", show_edges=True, cmap="rainbow", clim=[mn,mx])
    pl.add_mesh(beam, color="blue", line_width=5)
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
beamdisp   = np.loadtxt(f"{Dir}/BeamDisp0.out")
slicedmesh.points = originalcoords
beam.points = beamoriginalcoords

i = slicedmesh.find_closest_point([0,0,0.5])
i = beam.find_closest_point([0,0,3.0])

# plt.plot(time,array[:,i])
plt.plot(time,beamdisp[:,i*3])
plt.xlabel("Time [s]")
plt.ylabel("Displacement [m]")
plt.xlim([0,6])
plt.ylim([-0.02,0.02])
# plt.legend([f"{typename}"])
plt.title(f"Beam Displacement ({typename})")
plt.grid(linestyle="-.")
plt.show()



# %%
import numpy as np
import matplotlib.pyplot as plt 
x = np.linspace(0,140,280)



plt.plot(x,u)
# %%


typename = "PML"
# typename = "FIXED"
Dir = f"results/{typename}"
griddisp0   = np.loadtxt(f"{Dir}/NodeDisp0.out")
beamdisp0   = np.loadtxt(f"{Dir}/BeamDisp0.out")
gridaccel0  = np.loadtxt(f"{Dir}/NodeAccl0.out")
time0       = np.loadtxt(f"{Dir}/Time.out")[:,0]

typename = "FIXED"
# typename = "FIXED"
Dir = f"results/{typename}"
griddisp1   = np.loadtxt(f"{Dir}/NodeDisp0.out")
beamdisp1   = np.loadtxt(f"{Dir}/BeamDisp0.out")
gridaccel1  = np.loadtxt(f"{Dir}/NodeAccl0.out")
time1       = np.loadtxt(f"{Dir}/Time.out")[:,0]

# %%
i = beam.find_closest_point([0,0,3.0])
plt.plot(time1,beamdisp1[:,i*3],"b",label="Fixed boundary")
plt.plot(time0,beamdisp0[:,i*3],"r",label="PML boundary")
plt.xlabel("Time [s]")
plt.ylabel("Beam Head Displacement [m]")
plt.xlim([0,4])
plt.ylim([-0.02,0.02])
plt.legend(loc="upper right")
# %%
i = grid.find_closest_point([0,0,0.0])
plt.plot(time1,griddisp1[:,i*3],"b",label="Fixed boundary")
plt.plot(time0,griddisp0[:,i*3],"r",label="PML boundary")
plt.xlabel("Time [s]")
plt.ylabel("Foundation Displacement [m]")
plt.xlim([0,4])
plt.ylim([-0.02,0.02])
plt.legend(loc="upper right")
# %%
