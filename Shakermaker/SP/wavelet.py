# %%
import numpy as np
import matplotlib.pyplot as plt

def ricker(f, length=0.128, dt=0.001, center=0.0, amplitude=1.0):
    t = np.arange(-length/2, (length-dt)/2, dt)
    y = (1.0 - 2.0*(np.pi**2)*(f**2)*(t**2)) * np.exp(-(np.pi**2)*(f**2)*(t**2))
    t = t + length/2
    index2 = np.argmin(np.abs(t-length/2))
    index1 = np.argmin(np.abs(t-center))
    y = np.roll(y, -(index2 -index1))
    y = y * amplitude
    return t, y
 
f = 150 # A low wavelength of 25 Hz
t, w = ricker(f,1.0,0.0005,0.3, 25000)
plt.plot(t , w)
plt.title(f'Ricker Wavelet (Frequency: {f} Hz)')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.grid()

# save the force.dat
np.savetxt("force.dat", w, fmt="%.8f")

load = np.loadtxt("force.dat")
dt = 0.001
plt.plot(np.arange(0,0.0005*load.shape[0],0.0005), load)
plt.xlim(0.2, 0.6)
# %%
