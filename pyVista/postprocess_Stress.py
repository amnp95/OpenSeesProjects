# %%
import numpy as np 
import pandas as pd
import pyvista as pv
import os
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
# ==================================================================
# load nodes and elements dataframes
# ==================================================================
nodes = pd.read_hdf("nodes.h5")
elements = pd.read_hdf("elements.h5")     
pv.set_jupyter_backend('client')

# ==================================================================
# create PML mesh
# ==================================================================
nodes_in_mesh = nodes[nodes['Domain'] == "pml"]
elements_in_mesh = elements[elements['Domain'] == "pml"]
nodetoint = dict(zip(nodes_in_mesh['tag'], range(nodes_in_mesh.shape[0])))
cells = elements_in_mesh[['node1', 'node2', 'node3', 'node4', 'node5', 'node6', 'node7', 'node8']]
cells = cells.applymap(nodetoint.get)
cells = cells.to_numpy(dtype=int)
points = nodes_in_mesh[['x', 'y', 'z']].to_numpy(dtype=float)
celltypes = np.ones(cells.shape[0],dtype= int) * pv.CellType.HEXAHEDRON
cells = np.insert(cells, 0, 8, axis=1)
pmlgrid = pv.UnstructuredGrid(cells, celltypes.tolist(), points.tolist())
pmlgrid = pmlgrid.clip_box([-100,100,-100,100,-40, 0 ], invert=False,crinkle=True)
pmlgrid.clear_data()


# ==================================================================
# create whole grid mesh
# ==================================================================

nodes_in_mesh    =   nodes[nodes['Domain'] == "reg"]
elements_in_mesh = elements[elements['Domain'] == "reg"]
nodetoint = dict(zip(nodes_in_mesh['tag'], range(nodes_in_mesh.shape[0])))
cells = elements_in_mesh[['node1', 'node2', 'node3', 'node4', 'node5', 'node6', 'node7', 'node8']]
cells = cells.applymap(nodetoint.get)
cells = cells.to_numpy(dtype=int)
points = nodes_in_mesh[['x', 'y', 'z']].to_numpy(dtype=float)
celltypes = np.ones(cells.shape[0],dtype= int) * pv.CellType.HEXAHEDRON
cells = np.insert(cells, 0, 8, axis=1)
grid = pv.UnstructuredGrid(cells, celltypes.tolist(), points.tolist())


# ==================================================================
# create regular mesh
# ==================================================================
indexes  = grid.find_cells_within_bounds([-37.4,37.4,-1.25,1.3,-77.4,5.0])
reg      = grid.extract_cells(indexes)



# ==================================================================
# create DRM mesh
# ==================================================================
cube = pv.Cube(center=(0,0,0), x_length=80, y_length=80, z_length=55)
grid = grid.clip_box(cube,invert=False, crinkle=True)
cube = pv.Cube(center=(0,0,0), x_length=75, y_length=75, z_length=47.5)
drm = grid.clip_box(cube,invert=True, crinkle=True)
# indicies = drm.find_cells_within_bounds([-100,100,-1.25,1.3,-80,5.0])
# drm = drm.extract_cells(indicies)
reg = grid.clip_box(cube,invert=False, crinkle=True)

# ==================================================================
# create interface mesh
# ==================================================================
data = np.loadtxt("interfaceInfo.dat")
interface = pv.PolyData(data[:,1:4])



# create pile mesh 
# ==================================================================

zstart = -30
zend   = 20 
x      = 0.0
y      = 0.0
dz     = 0.5
pile   = pv.PolyData()
for z in np.arange(zstart,zend,dz):
    cube = pv.Cube(center=(x,y,z), x_length=0.5, y_length=0.5, z_length=dz)
    pile = pile + cube
pile = pile.clip_box([-100,100,-100,100,-20, 20 ], invert=False,crinkle=True)
pile.clear_data()


# create big mass add the head of the pile 
mass = pv.Cube(center=(0,0,zend+2.5), x_length=5, y_length=5, z_length=5)
mass = pv.Sphere(center=(0,0,zend+1.5), radius=2.5)

# set the camera position and width of the view window the same for all plots
# ==================================================================
pl = pv.Plotter(off_screen=False,window_size=[600,600])
pl.set_background('white')
# pl.set_position([0,0,0])
# pl.set_window_size(1920,1080)
# pl.set_scale(1,1,1)
# pl.set_focus([0,0,0])
# pl.camera_position = [(0, 0, 100), (0, 0, 0), (0, 1, 0)]
# pl.camera_set = True
# pl.camera = pl.camera
# pl.camera.zoom(1.5)
# matplotlib default color
# matplot lib defulat colors #1f77b4, #ff7f0e, #2ca02c, #d62728, #9467bd, #8c564b, #e377c2, #7f7f7f, #bcbd22, #17becf


pl.add_mesh(reg, show_edges=False, style="surface", opacity=0.7, color="#1f77b4")
pl.add_mesh(drm, show_edges=False, style="surface", opacity=1.0, color="#2ca02c")
pl.add_mesh(pmlgrid, show_edges=False, style="surface", opacity=1.0, color="#ff7f0e")
pl.add_mesh(interface, opacity=1.0, color="red", point_size=5)
pl.add_mesh(pile, show_edges=False, style="surface", opacity=1.0, color="blue")
pl.add_mesh(mass, show_edges=False, style="surface", opacity=1.0, color="blue")
# adding legend 
pl.add_legend(labels=[["Soil mesh", "#1f77b4"], ["DRM mesh", "#2ca02c"], ["PML mesh", "#ff7f0e"], ["Interface", "red"], ["Pile", "blue"]], border=True, bcolor="white", size=[0.3,0.3])
pl.show()



# # %%
# cube = pv.Cube(center=(0,0,30), x_length=20, y_length=20, z_length=60)
# # create a building mesh
# xrng = np.arange(-10, 10, 2.5, dtype=np.float32)
# yrng = np.arange(-10, 10, 2.5, dtype=np.float32)
# zrng = np.arange(0, 60, 5, dtype=np.float32)
# x, y, z = np.meshgrid(xrng, yrng, zrng, indexing='ij')
# cube = pv.StructuredGrid(x, y, z)

# zrange = np.arange(0, 2.1, 2., dtype=np.float32)
# x,y,z = np.meshgrid(xrng,yrng,zrange,indexing='ij')
# floor =pv.StructuredGrid(x,y,z)

# pl = pv.Plotter(off_screen=False)
# pl.set_background('white')
# # pl.add_mesh(cube, show_edges=True, style="wireframe", opacity=1.0, color="green", line_width=0.5)
# # pl.add_mesh(reg, show_edges=True, style="surface", opacity=1.0, color="#1f77b4")

# # make the points of the rego border red

# # crere columns of the building
# X , Y = np.meshgrid(xrng,yrng,indexing='ij')
# z = [0,55]
# for x in xrng[::2]:
#     for y in yrng[::2]:
#         clindeer = pv.Cylinder(center=(x,y,55/2.), direction=(0,0,1), radius=0.5, height=55)
#         # pl.add_mesh(pv.Line([x,y,z[0]],[x,y,z[1]]), color="blue", line_width=2.0)
#         pl.add_mesh(clindeer, show_edges=False, style="surface", opacity=1.0, color="red")

# for i in range(1,11):
#     f = floor.copy()
#     f.points[:,2] += 5  + 5 * i
#     pl.add_mesh(f, show_edges=False, style="surface", opacity=1.0)
    


# xrng = np.arange(-60, 40.1, 20, dtype=np.float32)
# yrng = np.arange(-60, 40.1, 20, dtype=np.float32)
# zrng = np.arange(-60, 0.1, 20, dtype=np.float32)
# x, y, z = np.meshgrid(xrng, yrng, zrng, indexing='ij')
# reg = pv.StructuredGrid(x, y, z)
# pl.add_mesh(reg, show_edges=False, style="Surface", opacity=0.25, color="purple",show_vertices=False)
# pl.add_mesh(reg.extract_feature_edges())
# # shrink the rego mesh
# xrng = np.arange(-30, 10.1, 10, dtype=np.float32)
# yrng = np.arange(-30, 10.1, 10, dtype=np.float32)
# zrng = np.arange(-30, 0.1, 10, dtype=np.float32)
# x, y, z = np.meshgrid(xrng, yrng, zrng, indexing='ij')
# reg2 = pv.StructuredGrid(x, y, z)
# pl.add_mesh(reg2, show_edges=False, style="Surface", opacity=0.5, color="#1f77b4",show_vertices=False,smooth_shading=True, split_sharp_edges=True)
# # show outer boundary
# surf = reg2.extract_surface()
# pl.add_mesh(surf, show_edges=False, style="points", opacity=1.0, color="red", point_size=5)
# pl.enable_ssao(kernel_size=512) 

# # adding legend
# pl.add_legend(labels=[["DRM Boundary", "#1f77b4"], ["HDF5 Container","pink"]], border=False, bcolor="white", size=[0.3,0.3],face="rectangle")


# # increase render quality
# pl.enable_eye_dome_lighting() 
# pl.renderer.use_depth_peeling=True
# pl.save_graphic("building.pdf")
# pl.screenshot("building.png")
# pl.show()

# #
# # %%
# #==================================================================
# from pyvista import examples
# structure, air = examples.download_electronics_cooling()
# # view plot
# jpeg = pv.read("build.png")
# # view the plot
# pl = pv.Plotter(off_screen=False)
# pl.add_mesh(reg, show_edges=True, style="surface", opacity=1.0, color="#1f77b4")
# pl.enable_ssao(kernel_size=128)
# pl.show()




# # %%

# %%
