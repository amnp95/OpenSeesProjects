#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft

def ricker_wavelet(t, f0, A=1):
    """
    Generate a Ricker wavelet pulse.

    Parameters:
        t (array): Time array.
        f0 (float): Central frequency of the wavelet.

    Returns:
        array: Ricker wavelet pulse.
    """
    A *= (1 - 2 * (np.pi * f0 * t)**2) * np.exp(-(np.pi * f0 * t)**2)
    return A


# Parameters
g     = 9.81    # gravitational acceleration
A     = 1*g    # Amplitude of the Ricker wavelet
dt    = 0.001   # Time step
t_max = 0.5     # Maximum time value
f0    = 10      # Central frequency of the Ricker wavelet
   

# Time array
t = np.arange(-t_max, t_max, dt)

# Generate Ricker wavelet pulse
wavelet = ricker_wavelet(t, f0, A)
t += 0.5


# plt.plot(t, wavelet)
# plt.xlabel("Time (s)")
# plt.ylabel("Amplitude")
# plt.xlim(0, 1.0)

# calculte velocity and displacement
velocity = np.cumsum(wavelet)*dt
displacement = np.cumsum(velocity)*dt


np.savetxt('ricker.acc', wavelet, fmt='%.7f')
np.savetxt('ricker.time', t, fmt='%.3f')
np.savetxt('ricker.vel', velocity, fmt='%.7f')
np.savetxt('ricker.dis', displacement, fmt='%.7f')


# %%
