# %%
import numpy as np
import matplotlib.pyplot as plt
import h5py
# %%
fig, ax = plt.subplots(3, 1, figsize=(10, 10),sharex=True)
data = np.loadtxt("NodeDisp.out")
controldispx = np.loadtxt("dispxcontrol.txt",delimiter=",")
controldispz = np.loadtxt("dispycontrol.txt",delimiter=",")

ax[0].plot(data[:, 0], data[:, [1,4,7]])
ax[0].plot(controldispx[0,:],controldispx[[1,2,3],:].T,"k--")
ax[1].plot(data[:, 0], data[:, [2,5,8]])
ax[2].plot(data[:, 0], data[:, [3,6,9]])
ax[2].plot(controldispz[0,:],controldispz[[1,2,3],:].T,"--")
ax[0].set_ylabel("X disp [m]")
ax[1].set_ylabel("Y disp [m]")
ax[2].set_ylabel("Z disp [m]")
ax[0].set_ylim([-0.011, 0.011])
ax[1].set_ylim([-0.011, 0.011])
ax[0].grid(True, which="both", ls="--")
ax[1].grid(True, which="both", ls="--")
ax[2].grid(True, which="both", ls="--")
ax[0].legend(["P1","P2","P3","analytical P1","analytical P2","analytical P3"], loc="upper right")
ax[2].legend(["P1","P2","P3","analytical P1","analytical P2","analytical P3"], loc="upper right")
fig.savefig("figs/fixed_Displacment_p1_p3.png")

# %%
fig, ax = plt.subplots(3, 1, figsize=(10, 10),sharex=True)
ax[0].plot(data[:, 0], data[:, [10,13,16]])
ax[0].plot(controldispx[0,:],controldispx[[4,5,6],:].T,"k--")
ax[1].plot(data[:, 0], data[:, [11,14,17]])
ax[2].plot(data[:, 0], data[:, [12,15,18]])
ax[2].plot(controldispz[0,:],controldispz[[4,5,6],:].T,"--")
ax[0].set_ylabel("X disp [m]")
ax[1].set_ylabel("Y disp [m]")
ax[2].set_ylabel("Z disp [m]")
# ax[0].set_ylim([-0.011, 0.011])
ax[1].set_ylim([-0.011, 0.011])
ax[0].grid(True, which="both", ls="--")
ax[1].grid(True, which="both", ls="--")
ax[2].grid(True, which="both", ls="--")
ax[0].legend(["P4","P5","P6","analytical P4","analytical P5","analytical P6"], loc="upper right")
ax[2].legend(["P4","P5","P6","analytical P4","analytical P5","analytical P6"], loc="upper right")
fig.savefig("figs/fixed_Displacment_p4_p6.png")
# %%
fig, ax = plt.subplots(3, 1, figsize=(10, 10),sharex=True)
data = np.loadtxt("NodeAccl.out")
controldispx = np.loadtxt("accxcontrol.txt",delimiter=",")
controldispz = np.loadtxt("accycontrol.txt",delimiter=",")

ax[0].plot(data[:, 0], data[:, [1,4,7]])
ax[0].plot(controldispx[0,:],controldispx[[1,2,3],:].T,"k--")
ax[1].plot(data[:, 0], data[:, [2,5,8]])
ax[2].plot(data[:, 0], data[:, [3,6,9]])
ax[2].plot(controldispz[0,:],controldispz[[1,2,3],:].T,"--")
ax[0].set_ylabel("X acc [m/s^2]")
ax[1].set_ylabel("Y acc [m/s^2]")
ax[2].set_ylabel("Z acc [m/s^2]")
# ax[0].set_ylim([-0.011, 0.011])
# ax[1].set_ylim([-0.011, 0.011])
ax[0].grid(True, which="both", ls="--")
ax[1].grid(True, which="both", ls="--")
ax[2].grid(True, which="both", ls="--")
ax[0].legend(["P1","P2","P3","analytical P1","analytical P2","analytical P3"], loc="upper right")
ax[2].legend(["P1","P2","P3","analytical P1","analytical P2","analytical P3"], loc="upper right")
fig.savefig("figs/fixed_Acceleration_p1_p3.png")
# %%
fig, ax = plt.subplots(3, 1, figsize=(10, 10),sharex=True)
ax[0].plot(data[:, 0], data[:, [10,13,16]])
ax[0].plot(controldispx[0,:],controldispx[[4,5,6],:].T,"k--")
ax[1].plot(data[:, 0], data[:, [11,14,17]])
ax[2].plot(data[:, 0], data[:, [12,15,18]])
ax[2].plot(controldispz[0,:],controldispz[[4,5,6],:].T,"--")
ax[0].set_ylabel("X acc [m/s^2]")
ax[1].set_ylabel("Y acc [m/s^2]")
ax[2].set_ylabel("Z acc [m/s^2]")
# ax[0].set_ylim([-0.011, 0.011])
ax[1].set_ylim([-0.011, 0.011])
ax[0].grid(True, which="both", ls="--")
ax[1].grid(True, which="both", ls="--")
ax[2].grid(True, which="both", ls="--")
ax[0].legend(["P4","P5","P6","analytical P4","analytical P5","analytical P6"], loc="upper right")
ax[2].legend(["P4","P5","P6","analytical P4","analytical P5","analytical P6"], loc="upper right")
fig.savefig("figs/fixed_Acceleration_p4_p6.png")

# %%
fig, ax = plt.subplots(3, 1, figsize=(10, 10),sharex=True)
data = np.loadtxt("NodeDisp_PML.out")
controldispx = np.loadtxt("dispxcontrol.txt",delimiter=",")
controldispz = np.loadtxt("dispycontrol.txt",delimiter=",")

ax[0].plot(data[:, 0], data[:, [1,4,7]])
ax[0].plot(controldispx[0,:],controldispx[[1,2,3],:].T,"k--")
ax[1].plot(data[:, 0], data[:, [2,5,8]])
ax[2].plot(data[:, 0], data[:, [3,6,9]])
ax[2].plot(controldispz[0,:],controldispz[[1,2,3],:].T,"--")
ax[0].set_ylabel("X disp [m]")
ax[1].set_ylabel("Y disp [m]")
ax[2].set_ylabel("Z disp [m]")
ax[0].set_ylim([-0.011, 0.011])
ax[1].set_ylim([-0.011, 0.011])
ax[0].grid(True, which="both", ls="--")
ax[1].grid(True, which="both", ls="--")
ax[2].grid(True, which="both", ls="--")
ax[0].legend(["P1","P2","P3","analytical P1","analytical P2","analytical P3"], loc="upper right")
ax[2].legend(["P1","P2","P3","analytical P1","analytical P2","analytical P3"], loc="upper right")
fig.savefig("figs/PML_Displacment_p1_p3.png")

# %%
fig, ax = plt.subplots(3, 1, figsize=(10, 10),sharex=True)
ax[0].plot(data[:, 0], data[:, [10,13,16]])
ax[0].plot(controldispx[0,:],controldispx[[4,5,6],:].T,"k--")
ax[1].plot(data[:, 0], data[:, [11,14,17]])
ax[2].plot(data[:, 0], data[:, [12,15,18]])
ax[2].plot(controldispz[0,:],controldispz[[4,5,6],:].T,"--")
ax[0].set_ylabel("X disp [m]")
ax[1].set_ylabel("Y disp [m]")
ax[2].set_ylabel("Z disp [m]")
# ax[0].set_ylim([-0.011, 0.011])
ax[1].set_ylim([-0.011, 0.011])
ax[0].grid(True, which="both", ls="--")
ax[1].grid(True, which="both", ls="--")
ax[2].grid(True, which="both", ls="--")
ax[0].legend(["P4","P5","P6","analytical P4","analytical P5","analytical P6"], loc="upper right")
ax[2].legend(["P4","P5","P6","analytical P4","analytical P5","analytical P6"], loc="upper right")
fig.savefig("figs/PML_Displacment_p4_p6.png")
# %%
fig, ax = plt.subplots(3, 1, figsize=(10, 10),sharex=True)
data = np.loadtxt("NodeAccl_PML.out")
controldispx = np.loadtxt("accxcontrol.txt",delimiter=",")
controldispz = np.loadtxt("accycontrol.txt",delimiter=",")

ax[0].plot(data[:, 0], data[:, [1,4,7]])
ax[0].plot(controldispx[0,:],controldispx[[1,2,3],:].T,"k--")
ax[1].plot(data[:, 0], data[:, [2,5,8]])
ax[2].plot(data[:, 0], data[:, [3,6,9]])
ax[2].plot(controldispz[0,:],controldispz[[1,2,3],:].T,"--")
ax[0].set_ylabel("X acc [m/s^2]")
ax[1].set_ylabel("Y acc [m/s^2]")
ax[2].set_ylabel("Z acc [m/s^2]")
# ax[0].set_ylim([-0.011, 0.011])
# ax[1].set_ylim([-0.011, 0.011])
ax[0].grid(True, which="both", ls="--")
ax[1].grid(True, which="both", ls="--")
ax[2].grid(True, which="both", ls="--")
ax[0].legend(["P1","P2","P3","analytical P1","analytical P2","analytical P3"], loc="upper right")
ax[2].legend(["P1","P2","P3","analytical P1","analytical P2","analytical P3"], loc="upper right")
fig.savefig("figs/PML_Acceleration_p1_p3.png")
# %%
fig, ax = plt.subplots(3, 1, figsize=(10, 10),sharex=True)
ax[0].plot(data[:, 0], data[:, [10,13,16]])
ax[0].plot(controldispx[0,:],controldispx[[4,5,6],:].T,"k--")
ax[1].plot(data[:, 0], data[:, [11,14,17]])
ax[2].plot(data[:, 0], data[:, [12,15,18]])
ax[2].plot(controldispz[0,:],controldispz[[4,5,6],:].T,"--")
ax[0].set_ylabel("X acc [m/s^2]")
ax[1].set_ylabel("Y acc [m/s^2]")
ax[2].set_ylabel("Z acc [m/s^2]")
# ax[0].set_ylim([-0.011, 0.011])
ax[1].set_ylim([-0.011, 0.011])
ax[0].grid(True, which="both", ls="--")
ax[1].grid(True, which="both", ls="--")
ax[2].grid(True, which="both", ls="--")
ax[0].legend(["P4","P5","P6","analytical P4","analytical P5","analytical P6"], loc="upper right")
ax[2].legend(["P4","P5","P6","analytical P4","analytical P5","analytical P6"], loc="upper right")
fig.savefig("figs/PML_Acceleration_p4_p6.png")
# %%
