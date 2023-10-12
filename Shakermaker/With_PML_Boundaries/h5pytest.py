# %%
import h5py
import numpy as np
import matplotlib.pyplot as plt

# f = h5py.File('test.h5drm', 'r')
f = h5py.File(name='/home/amin/Projects/Github/OpenSeesProjects/Shakermaker/With_PML_Boundaries/test.h5drm', mode='r')
acc_center = np.array(f["DRM_QA_Data"]["acceleration"])

f["DRM_Metadata"].keys()
dt = np.array(f["DRM_Metadata"]["dt"]).__float__()
tend = np.array(f["DRM_Metadata"]["tend"]).__float__()
tstart = np.array(f["DRM_Metadata"]["tstart"]).__float__()


# %%
fig,ax = plt.subplots(3,1,sharex=True,sharey=True)
time = np.arange(tstart,tend,dt)
ax[0].plot(time ,acc_center[0,:],label="Shaker maker")
ax[1].plot(time ,acc_center[1,:])
ax[2].plot(time ,acc_center[2,:])
ax[2].set_xlim([10,20])

# data = np.loadtxt("/home/amin/Projects/Github/OpenSeesProjects/Shakermaker/With_PML_Boundaries/NodeAccl0.out")
# time = data[:,0]
# ax[0].plot(time ,data[:,1],label= "abaqus")
# ax[1].plot(time ,data[:,2])
# ax[2].plot(time ,data[:,3])
ax[0].legend()
ax[0].set_ylabel('Acceleration x [m/s]')
ax[1].set_ylabel('Acceleration y [m/s]')
ax[2].set_ylabel('Acceleration z [m/s]')
# ax[0].set_ylim(-0.0175,0.01)
plt.show()
# %%

# %%
