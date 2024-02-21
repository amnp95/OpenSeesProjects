#%%
import numpy as np
import pyvista as pv
import time
#%%
time1 = time.time()
# create a strucrured grid 
# xwidth  = 500
# ywidth  = 500
# zwidth  = 500
# eps     = 1e-6

# Create Davis Centrifuge container
xwidth  = 0.6
ywidth  = 1.7
zwidth  = 0.7
eps     = 1e-6

# Xmeshsize, Ymeshsize, Zmeshsize = (2.5, 2.5, 2.5)
# Xmeshsize, Ymeshsize, Zmeshsize = (5, 5, 5)

# Davis Centrifuge container
Xmeshsize, Ymeshsize, Zmeshsize = (0.025, 0.025, 0.025)

x = np.arange(-xwidth/2., xwidth/2.+eps, Xmeshsize)
y = np.arange(-ywidth/2., ywidth/2.+eps, Ymeshsize)
z = np.arange(-zwidth, 0+eps, Zmeshsize)

x, y, z = np.meshgrid(x, y, z)

mesh = pv.StructuredGrid(x, y, z)


PMLthickness = 5.0
PMLXmeshsize, PMLYmeshsize, PMLZmeshsize = (2.5, 2.5, 2.5)

# Davis centrifuge containes
PMLthickness = 0.07
PMLXmeshsize, PMLYmeshsize, PMLZmeshsize = (0.035, 0.035, 0.035)

# sperate PML layer 
xmin = -xwidth/2. + PMLthickness
xmax = xwidth/2. - PMLthickness
ymin = -ywidth/2. + PMLthickness
ymax = ywidth/2. - PMLthickness
zmin = -zwidth + PMLthickness
zmax = 0
cube = pv.Cube(bounds=[xmin,xmax,ymin,ymax,zmin,zmax])
PML = mesh.clip_box(cube,invert=True,crinkle=True,progress_bar = True)
PML.cell_data["PML"] = np.ones(PML.n_cells,dtype=bool)


# seperate regular mesh
reg = mesh.clip_box(cube,invert=False,crinkle=True,progress_bar = True)
reg.cell_data["PML"] = np.zeros(reg.n_cells,dtype=bool)

numreg = reg.n_cells
numpml = PML.n_cells

mesh = reg.merge(PML,merge_points=False,tolerance=1e-6,progress_bar = True)
print(mesh.n_cells,mesh.n_points)
mapping = mesh.clean(produce_merge_map=True)["PointMergeMap"]

regindicies = np.where(mapping[PML.n_points:]<PML.n_points)[0]
PMLindicies = mapping[PML.n_points+regindicies]

mesh.point_data["boundary"] = np.zeros(mesh.n_points,dtype=int)-1
mesh.point_data["boundary"][PMLindicies] = regindicies + PML.n_points
mesh.point_data["boundary"][PML.n_points + regindicies] = PMLindicies 

time2 = time.time()
print(time2-time1)


# %%

mesh["boundary"] = np.where(mesh.point_data["boundary"]>0,True,False)
pl = pv.Plotter()
pl.add_mesh(mesh,scalars="boundary",show_edges=True, style="surface", opacity=1.0, color="#2ca02c")
pl.show()


# pl.add_mesh(drm, show_edges=False, style="surface", opacity=1.0, color="#2ca02c")



# %%
