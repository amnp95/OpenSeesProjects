# %%
import numpy as np
import matplotlib.pyplot as plt
import h5py
data = np.loadtxt('NodeAccl1.out')
fig, ax = plt.subplots(3, 1, figsize=(8, 8),sharex=True, sharey=True)
ax[0].plot(data[:, 0], data[:, 1], label="OpenSees")
ax[1].plot(data[:, 0], data[:, 2])
ax[2].plot(data[:, 0], data[:, 3])
ax[0].set_ylabel('Acceleration x [m/s]')
ax[1].set_ylabel('Acceleration y [m/s]')
ax[2].set_ylabel('Acceleration z [m/s]')
ax[0].grid(linestyle='--')
ax[1].grid(linestyle='--')
ax[2].grid(linestyle='--')
ax[2].set_xlim([8, 20])
# ax[2].set_ylim([-4e10, 5e10])
plt.xlabel('Time (s)')
# plt.ylim(-0.06, 0.04)
# plt.xlim(right=1.1)
# centacc = np.loadtxt("acceleration.acc")
# ax[2].plot(centacc[:, 0], centacc[:, 1]*5)
# ax[0].plot(centacc[:, 0], centacc[:, 2]*5, label='Shakermaker')
# ax[1].plot(centacc[:, 0], centacc[:, 3]*5)
f = h5py.File(name='/home/amin/Projects/Github/OpenSeesProjects/Shakermaker/With_PML_Boundaries/test.h5drm', mode='r')
acc_center = np.array(f["DRM_QA_Data"]["acceleration"])

f["DRM_Metadata"].keys()
dt = np.array(f["DRM_Metadata"]["dt"]).__float__()
tend = np.array(f["DRM_Metadata"]["tend"]).__float__()
tstart = np.array(f["DRM_Metadata"]["tstart"]).__float__()
time = np.arange(tstart,tend,dt)
ax[0].plot(time ,acc_center[0,:],label="Shaker maker")
ax[1].plot(time ,acc_center[1,:])
ax[2].plot(time ,acc_center[2,:])
ax[0].legend()
plt.show()

# %%
# data = np.loadtxt('NodeAccl0.out')
# time = data[:, 0]
# acc = data[:, 1:]
# fig, ax = plt.subplots(3, 1, figsize=(8, 8),sharex=True, sharey=True)

# ax[0].plot(time, acc[:, 0], label='x')
# ax[1].plot(time, acc[:, 1], label='y')
# ax[2].plot(time, acc[:, 2], label='z')

# ax[0].set_ylabel('Acceleration x [m/s^2]')
# ax[1].set_ylabel('Acceleration y [m/s^2]')
# ax[2].set_ylabel('Acceleration z [m/s^2]')
# ax[2].set_xlim([13, 20])
# ax[2].set_xlabel('Time [s]')


# centacc = np.loadtxt("acceleration.acc")
# ax[2].plot(centacc[:, 0], centacc[:, 1], label='x')
# ax[0].plot(centacc[:, 0], centacc[:, 2], label='y')
# ax[1].plot(centacc[:, 0], centacc[:, 3], label='z')

# %%
