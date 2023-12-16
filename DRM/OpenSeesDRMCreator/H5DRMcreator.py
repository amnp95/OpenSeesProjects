# %%
import numpy as np
import h5py

filename     = "DRMload.h5drm"
coords       = np.loadtxt("coords.txt",       comments="#" , delimiter=",")
displacement = np.loadtxt("disp.txt", comments="#" , delimiter=",")
acceleration = np.loadtxt("acc.txt", comments="#" , delimiter=",")
velocity     = np.loadtxt("vel.txt",     comments="#" , delimiter=",")
nstations    = coords.shape[0]
spacing      = [5, 5, 5]
Dimensions   = [(25,-25), (25,25), (0,20)]
dt           = 0.01
tsart        = 0
name         = "test"
# %%
# loop over the corrds to say that they are internal or not
xmax, xmin = Dimensions[0]
ymax, ymin = Dimensions[1]
zmax, zmin = Dimensions[2]

internal = np.ones(nstations, dtype=np.bool)
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

# %%
# read hdf5 file
HDffile = h5py.File("test.h5drm", mode="r")
HDffile.keys()
HDffile["/DRM_Metadata"].keys()
HDffile["/DRM_Metadata"]["drmbox_xmin"][()]
HDffile["/DRM_Metadata"]["drmbox_xmax"][()]
HDffile["/DRM_Metadata"]["drmbox_ymin"][()]
HDffile["/DRM_Metadata"]["drmbox_ymax"][()]
HDffile["/DRM_Metadata"]["drmbox_zmin"][()]
HDffile["/DRM_Metadata"]["drmbox_zmax"][()]
HDffile["/DRM_Metadata"]["dt"][()]
HDffile["/DRM_Metadata"]["tend"][()]
HDffile["/DRM_Metadata"]["tstart"][()]
HDffile["/DRM_Metadata"]["program_used"][()]
HDffile["/DRM_Metadata"]["h"][()]
# HDffile["/DRM_Metadata"]["name"][()]

# %%
# open a file in writing mode
DRMfile = h5py.File(filename, mode="w")


# create gropus
DRMdata = DRMfile.create_group("/DRM_Data")
DRMfile.create_group("//DRM_QA_Data")
DRMmetadata = DRMfile.create_group("/DRM_Metadata")


# creating dataset
DRMdata.create_dataset("xyz",          data=coords.T,       dtype=np.double)
DRMdata.create_dataset("internal",     data=internal,       dtype=bool)
DRMdata.create_dataset("displacement", data=displacement.T, dtype=np.double)
DRMdata.create_dataset("acceleration", data=acceleration.T, dtype=np.double)
DRMdata.create_dataset("velocity",     data=velocity.T,     dtype=np.double)
DRMdata.create_dataset("data_location",data=data_location)

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

# close the file
DRMfile.close()





# %%
