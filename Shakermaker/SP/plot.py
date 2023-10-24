# %%
import numpy as np
import matplotlib.pyplot as plt
import h5py
# %% 
# read hpy file
f = h5py.File("test.h5drm", "r")
acc = np.array(f["DRM_QA_Data"]["acceleration"])

# dt = np.array(f["DRM_Metadata"]["dt"])

tend = np.array(f["DRM_Metadata"]["tend"]).item()
tstart = np.array(f["DRM_Metadata"]["tstart"]).item()
dt = np.array(f["DRM_Metadata"]["dt"]).item()
time = np.arange(tstart, tend, dt)
acc.shape
# %%
fig, ax = plt.subplots(3, 1, figsize=(10, 10),sharex=True, sharey=True)

ax[0].plot(time,acc[0,:])
ax[1].plot(time,acc[1,:])
ax[2].plot(time,acc[2,:], label="Shaker maker")


# data = np.loadtxt("NodeAccl.out")
# ax[0].plot(data[:, 0], data[:, 1])
# ax[1].plot(data[:, 0], data[:, 2])
# ax[2].plot(data[:, 0], data[:, 3], label="Fixed boundaries")

# data = np.loadtxt("NodeAccl_damping.out")
# ax[0].plot(data[:, 0], data[:, 1])
# ax[1].plot(data[:, 0], data[:, 2])
# ax[2].plot(data[:, 0], data[:, 3], label="Numerical damping")


# data = np.loadtxt("NodeAccl_reyleigh.out")
# ax[0].plot(data[:, 0], data[:, 1])
# ax[1].plot(data[:, 0], data[:, 2])
# ax[2].plot(data[:, 0], data[:, 3], label="reyleigh damping")


# data = np.loadtxt("NodeAccl_PML1.out")
# ax[0].plot(data[:, 0], data[:, 1])
# ax[1].plot(data[:, 0], data[:, 2])
# ax[2].plot(data[:, 0], data[:, 3], label="PML")


data = np.loadtxt("NodeAccl_PML.out")
ax[0].plot(data[:, 0], data[:, 1])
ax[1].plot(data[:, 0], data[:, 2])
ax[2].plot(data[:, 0], data[:, 3], label="PML")
# data = np.loadtxt("NodeAccl_PML_M.out")
# ax[0].plot(data[:, 0], data[:, 1])
# ax[1].plot(data[:, 0], data[:, 2])
# ax[2].plot(data[:, 0], data[:, 3], label="without M")
# data = np.loadtxt("NodeAccl_PML_G.out")
# ax[0].plot(data[:, 0], data[:, 1])
# ax[1].plot(data[:, 0], data[:, 2])
# ax[2].plot(data[:, 0], data[:, 3], label="without G")

ax[2].set_xlim(8, 25)



# ax[2].set_ylim(-0.06, 0.05)


ax[2].legend(loc=1)




# %%
