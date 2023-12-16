# %%
import numpy as np
import h5py
from datetime import date

filename     = "SurfaceWave.h5drm"
author       = "Amin Pakzad"
daate        =  date.today()
coords       = np.loadtxt("coords.txt", comments="#" , delimiter=",")
displacement = np.loadtxt("disp.txt",   comments="#" , delimiter=",")
acceleration = np.loadtxt("acc.txt",    comments="#" , delimiter=",")
velocity     = np.loadtxt("vel.txt",    comments="#" , delimiter=",")
#%%
nstations    = coords.shape[0]
x0           = [0., 0., 0.]
spacing      = [0.5, 0.5, 0.5]
Dimensions   = [(9.5,-9.5), (9.5,-9.5), (0,-9.5)]
dt           = 0.0001
tsart        = 0.
name         = "SurfaceWave"
# %%
# scaling the time step 
scale = 1
if scale > 1:
    dt = dt * scale
    displacement = displacement[:,::scale]
    acceleration = acceleration[:,::scale]
    velocity     = velocity[:,::scale]
# %%
# loop over the corrds to say that they are internal or not
xmax, xmin = Dimensions[0]
ymax, ymin = Dimensions[1]
zmax, zmin = Dimensions[2]

internal = np.ones(nstations, dtype=bool)
tol = 1e-6
for i in range(nstations):
    x, y, z = coords[i, :]
    if x < xmin - tol or x > xmax + tol:
        internal[i] = False
        continue
    if y < ymin - tol or y > ymax + tol:
        internal[i] = False
        continue
    if z < zmin - tol or z > zmax + tol:
        internal[i] = False
        continue

data_location = np.arange(0, nstations, dtype=np.int32) * 3
xmax = xmax + spacing[0]
xmin = xmin - spacing[0]
ymax = ymax + spacing[1]
ymin = ymin - spacing[1]
zmax = zmax
zmin = zmin - spacing[2]
tend = dt * (displacement.shape[1] - 1) + tsart
program_used = "OpenSees"

# # %%
# # plot the nodes in 3d and the internal nodes with different color 
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# idx = np.where(internal == True)[0]
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')  # Use add_subplot instead of Axes3D
# ax.scatter(coords[internal, 0], coords[internal, 1], coords[internal, 2],color='b',s = 0.5)
# ax.scatter(coords[~internal, 0], coords[~internal, 1], coords[~internal, 2],color='r',s = 0.1)
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')
# plt.show()

# %%
# open a file in writing mode
DRMfile = h5py.File(filename, mode="w")


# create gropus
DRMdata = DRMfile.create_group("/DRM_Data")
DRMfile.create_group("//DRM_QA_Data")
DRMmetadata = DRMfile.create_group("/DRM_Metadata")


# creating dataset
DRMdata.create_dataset("xyz",          data=coords,         dtype=np.double)
DRMdata.create_dataset("internal",     data=internal,       dtype=bool)
DRMdata.create_dataset("displacement", data=displacement, dtype=np.double)
DRMdata.create_dataset("acceleration", data=acceleration, dtype=np.double)
DRMdata.create_dataset("velocity",     data=velocity,     dtype=np.double)
DRMdata.create_dataset("data_location",data=data_location)

DRMmetadata.create_dataset("drmbox_x0",   data=x0)
DRMmetadata.create_dataset("drmbox_xmin", data=xmin)
DRMmetadata.create_dataset("drmbox_xmax", data=xmax)
DRMmetadata.create_dataset("drmbox_ymin", data=ymin)
DRMmetadata.create_dataset("drmbox_ymax", data=ymax)
DRMmetadata.create_dataset("drmbox_zmin", data=zmin)
DRMmetadata.create_dataset("drmbox_zmax", data=zmax)
DRMmetadata.create_dataset("dt", data=dt)
DRMmetadata.create_dataset("tend", data=tend)
DRMmetadata.create_dataset("tstart", data=tsart)
DRMmetadata.create_dataset("program_used", data=program_used)
DRMmetadata.create_dataset("h", data=spacing)
DRMmetadata.create_dataset("name", data=name)
# DRMmetadata.create_dataset("created_by", data=author)
# DRMmetadata.create_dataset("created_on", data=daate)


# close the file
DRMfile.close()





# %%
