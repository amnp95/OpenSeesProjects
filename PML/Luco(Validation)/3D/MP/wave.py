# %%
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
A     = 857.32  # Amplitude of the Ricker wavelet
dt    = 0.001   # Time step
t_max = 0.5     # Maximum time value
f0    = 8      # Central frequency of the Ricker wavelet
# f0    = 2.83 
# omegaL = 14.69
# omegaH = 36.74
# f0 = omegaL/(2*np.pi)
# f0 = omegaH/(2*np.pi)   

# Time array
t = np.arange(-t_max, t_max, dt)

# Generate Ricker wavelet pulse
wavelet = ricker_wavelet(t, f0, A)
t += 0.5


# save the wavelet in file force.dat with float format
np.savetxt('force.dat', wavelet, fmt='%.6f')


# Time-domain plot
plt.figure(figsize=(10, 5))
plt.plot(t, wavelet, label='Ricker Wavelet')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Ricker Wavelet Pulse (Time Domain)')
# data = np.loadtxt("RichterHF.txt")
# plt.plot(np.arange(0, 1, 0.001), data, label='RichterLF')
plt.legend()
plt.grid(True)
plt.show()

# Frequency-domain plot
N = len(wavelet)
freq = np.fft.fftfreq(N, dt)
fft_wavelet = fft(wavelet)

plt.figure(figsize=(10, 5))
plt.plot(freq[0:100], np.abs(fft_wavelet)[0:100],"o--", label='Frequency Spectrum')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.title('Ricker Wavelet Pulse (Frequency Domain)')
plt.legend()
plt.grid(True)
plt.xlim(0, 2 * f0+10)  # Plot only positive frequencies up to 2 times the central frequency
plt.show()

# %%
