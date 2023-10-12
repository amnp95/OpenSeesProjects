from matplotlib.pylab import *
import h5py
d_rec = loadtxt("./analysis/results-dis-1.txt")
a_rec = loadtxt("./analysis/results-acc-1.txt")
Nt_rec = d_rec.shape[0]
dt_rec = 0.001
t_rec = arange(Nt_rec)*dt_rec


h5 = h5py.File("test.h5drm")
d_drm = h5["/DRM_QA_Data/displacement"][:].T
a_drm = h5["/DRM_QA_Data/acceleration"][:].T
Nt_drm  = d_drm.shape[0]
dt_drm = 0.01
t_drm = arange(Nt_drm)*dt_drm

print(f"{d_drm=}")
print(f"{Nt_drm=}")

figure(1)

ax0 = subplot(3,1,1)
subplot(3,1,2,sharex=ax0,sharey=ax0)
subplot(3,1,3,sharex=ax0,sharey=ax0)


for comp in range(3):
    subplot(3,1,1+comp)
    plot(t_rec, d_rec[:,comp],label="OpenSees")
    plot(t_drm, d_drm[:,comp],label="ShakerMaker")
    ylabel(f"Dis {['X','Y','Z'][comp]}")
legend()
xlabel("Time (s)")


figure(2)
ax0 = subplot(3,1,1)
subplot(3,1,2,sharex=ax0,sharey=ax0)
subplot(3,1,3,sharex=ax0,sharey=ax0)
for comp in range(3):
    subplot(3,1,1+comp)
    plot(t_rec, a_rec[:,comp],label="OpenSees")
    plot(t_drm, a_drm[:,comp],label="ShakerMaker")
    ylabel(f"Acc {['X','Y','Z'][comp]}")
legend()
xlabel("Time (s)")


show()