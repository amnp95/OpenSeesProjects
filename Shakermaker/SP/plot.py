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
# %%
fig, ax = plt.subplots(3, 1, figsize=(8, 8),sharex=True, sharey=True)

data = np.loadtxt("NodeAccl.out")
ax[0].plot(data[:, 0], data[:, 1], label="Fixed boundaries")
ax[1].plot(data[:, 0], data[:, 2])
ax[2].plot(data[:, 0], data[:, 3])
ax[0].set_xlim(0, 1.0)
data = np.loadtxt("NodeAccl_PML.out")
ax[0].plot(data[:, 0], data[:, 1], label="PML boundaries")
ax[1].plot(data[:, 0], data[:, 2])
ax[2].plot(data[:, 0], data[:, 3])
ax[0].legend()
# %%
# calcute the frequency content of the displacement
disp = data[:,3]
dt = 0.001
n = 1024
n = disp.shape[0]
disp = disp[:n]
freq = np.fft.fftfreq(n, dt)[:int(n/2)]
fft = np.fft.fft(disp)[:int(n/2)]
plt.plot(freq, np.abs(fft),"o-", label="FFT of disp")
plt.xlim(0,20)
plt.grid(linestyle="--")
# %%


# data = np.loadtxt("NodeAccl_PML.out")
# ax[0].plot(data[:, 0], data[:, 1], label="PML boundaries")
# ax[1].plot(data[:, 0], data[:, 2])
# ax[2].plot(data[:, 0], data[:, 3])
# ax[0].plot(time,acc[0,:], label="Shaker maker")
# ax[1].plot(time,acc[1,:])

# # ax[2].plot(time,acc[2,:])
# file = open("5566AXYZ_Abaqus.txt", "r")

# # read all lines
# lines = file.readlines()


# data = []
# for line in lines :
#     # delete the spaces
#     line = line.strip()
#     data.append(line.split())

# data = np.array(data, dtype=float)
# ax[0].plot(data[:,0], data[:,1], label="Abaqus")
# ax[1].plot(data[:,0], data[:,2])
# ax[2].plot(data[:,0], -data[:,3])

# file = open("5566AXYZ_Abaqus_1m .rpt", "r")

# # read all lines
# lines = file.readlines()


# data = []
# for line in lines :
#     # delete the spaces
#     line = line.strip()
#     data.append(line.split())

# data = np.array(data, dtype=float)
# ax[0].plot(data[:,0], data[:,1], label="Abaqus_1m")
# ax[1].plot(data[:,0], data[:,2])
# ax[2].plot(data[:,0], -data[:,3])


ax[0].set_xlim(0, 1.2)
ax[2].legend()
    # %%
fig, ax = plt.subplots()


data = np.loadtxt("force.dat")
# ax.plot(np.arange(0,0.001*data.shape[0],0.001), -data, label="Force")

# plot the frequncy content of the force 
force = data
dt = 0.001
n = force.shape[0]
freq = np.fft.fftfreq(n, dt)
fft = np.fft.fft(force)
ax.plot(freq, np.abs(fft), label="FFT of force")
ax.set_xlim(0, 12)



# %%
file = open("5566AXYZ_Abaqus.txt", "r")

# read all lines
lines = file.readlines()


data = []
for line in lines :
    # delete the spaces
    line = line.strip()
    data.append(line.split())

data = np.array(data, dtype=float)
data


# %%

VS = 63
H = 20

for n in range(0,10):
    w= VS/H * (3.1415926/2 + n*3.1415926)
    print(f"n={n}, w={w}, f={w/(2*3.1415926)}")
# %%