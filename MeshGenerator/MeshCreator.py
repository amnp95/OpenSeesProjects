#%%
import os 
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pyvista as pv
from Partition import partition
import matplotlib.pyplot as plt
import time

# chenge the directory to the current file



# create a dictionary for meshing information
info = {
    "RegularDomain": 1,
    "DRMDomain": 2,
    "PMLDomain": 3,
}

time1 = time.time()
# create a strucrured grid 
xwidth  = 20
ywidth  = 20
zwidth  = 50
eps     = 1e-6

Xmeshsize, Ymeshsize, Zmeshsize = (1, 1, 1)

x = np.arange(-xwidth/2., xwidth/2.+eps, Xmeshsize)
y = np.arange(-ywidth/2., ywidth/2.+eps, Ymeshsize)
z = np.arange(-zwidth, 0+eps, Zmeshsize)

x, y, z = np.meshgrid(x, y, z)

mesh = pv.StructuredGrid(x, y, z)


PMLthickness = 2.
PMLXmeshsize, PMLYmeshsize, PMLZmeshsize = (1., 1., 1.)

# sperate PML layer 
xmin = -xwidth/2. + PMLthickness
xmax = xwidth/2. - PMLthickness
ymin = -ywidth/2. + PMLthickness
ymax = ywidth/2. - PMLthickness
zmin = -zwidth + PMLthickness
zmax = 0
cube = pv.Cube(bounds=[xmin,xmax,ymin,ymax,zmin,zmax])
PML = mesh.clip_box(cube,invert=True,crinkle=True,progress_bar = True)
reg = mesh.clip_box(cube,invert=False,crinkle=True,progress_bar = True)


# now find DRM layer
indices = reg.find_cells_within_bounds([xmin + Xmeshsize + eps,
                              xmax - Xmeshsize - eps,
                              ymin + Ymeshsize + eps,
                              ymax - Ymeshsize - eps,
                              zmin + Zmeshsize + eps,
                              zmax + Zmeshsize + eps])

# now create complemntary indices for DRM
DRMindices = np.ones(reg.n_cells,dtype=bool)
DRMindices[indices] = False
DRMindices = np.where(DRMindices)[0]



reg.cell_data['Domain'] = np.ones(reg.n_cells,dtype=np.int8)*info["DRMDomain"]
reg.cell_data['Domain'][indices] = info["RegularDomain"]
PML.cell_data['Domain'] = np.ones(PML.n_cells,dtype=np.int8)*info["PMLDomain"]
reg.cell_data['partitioned'] = np.zeros(reg.n_cells,dtype=np.int32)


# partitioning regular mesh
regular = reg.extract_cells(indices,progress_bar=True)


DRM = reg.extract_cells(DRMindices,progress_bar=True)

reg_num_cores = 4
DRM_num_cores = 2
PML_num_cores = 3



partition(regular,reg_num_cores)
partition(DRM,DRM_num_cores)
partition(PML,PML_num_cores)

reg.cell_data['partitioned'][regular["vtkOriginalCellIds"]] = regular.cell_data['partitioned']
reg.cell_data['partitioned'][DRM["vtkOriginalCellIds"]] = DRM.cell_data['partitioned'] + reg_num_cores
PML.cell_data['partitioned'] = PML.cell_data['partitioned'] + reg_num_cores + DRM_num_cores


# merging PML and regular mesh to create a single mesh
mesh = reg.merge(PML,merge_points=False,tolerance=1e-6,progress_bar = True)


# mapping between PML and regular mesh on the boundary
mapping = mesh.clean(produce_merge_map=True)["PointMergeMap"]
regindicies = np.where(mapping[PML.n_points:]<PML.n_points)[0]
PMLindicies = mapping[PML.n_points+regindicies]


mesh.point_data["boundary"] = np.zeros(mesh.n_points,dtype=int)-1
mesh.point_data["boundary"][PMLindicies] = regindicies + PML.n_points
mesh.point_data["boundary"][PML.n_points + regindicies] = PMLindicies 

indices = np.where(mesh.point_data["boundary"]>0)[0]

time2 = time.time()
print("Time for meshing: ",time2-time1)

pl = pv.Plotter()
matplotlib_defaultcolors = [format(c) for c in plt.rcParams['axes.prop_cycle'].by_key()['color']]
pl.add_mesh(mesh, scalars="partitioned", show_edges=True,opacity=1.0,cmap=matplotlib_defaultcolors)
pl.show()


# download a topograpy data from a base map area of 10*10 km 
# https://portal.opentopography.org/raster?opentopoID=OTDS.082020.32611.1


num_layers = 3








# # %%
DRM.plot(scalars="partitioned", show_edges=True,opacity=1.0,cmap=matplotlib_defaultcolors)
# # %%
# PML.plot(scalars="partitioned", show_edges=True,opacity=1.0,cmap=matplotlib_defaultcolors)
# # %%
# regular.plot(scalars="partitioned", show_edges=True,opacity=1.0,cmap=matplotlib_defaultcolors)

# %%
# now writung the mesh to a file
# create the 















# %%

Dir = "OpenSeesMesh"
if not os.path.exists(Dir):
    os.makedirs(Dir)


min_core = mesh.cell_data['partitioned'].min()
max_core = mesh.cell_data['partitioned'].max()

# write the  mesh nodes
for core  in range(min_core,max_core+1):
    tmp  = mesh.extract_cells(np.where(mesh.cell_data['partitioned']==core)[0])
    f  = open(Dir + "/Nodes_" + str(core) + ".tcl", "w")

    for i in range(tmp.n_points):
        f.write(f"node  {tmp['vtkOriginalPointIds'][i]} {tmp.points[i][0]} {tmp.points[i][1]} {tmp.points[i][2]}\n")
    f.close()


# writing the mesh elements
f  = open(Dir + "/Elements.tcl", "w")
for eletag in range(mesh.n_cells):
    f.write(f"element {eletag}  {' '.join(str(x) for x in mesh.get_cell(eletag).point_ids)}\n")
f.close()












# %%
