from matplotlib.pylab import *

Vs = 200.
nu = 0.25
rho = 200.
fmax = 20
tmax = 20.


G = rho * Vs**2
E = G*(2*(1+nu))
M = 2*G*(1-nu)/(1-2*nu)
Vp = sqrt(M/rho)



print(f"{Vs=} {nu=} {rho=} {G=} {E=} {M=} {Vp=}")



λ = Vs/fmax

dx = λ/20

dt = dx / Vp
dt = 0.001

print(f"{dt = } < {dx / Vp=} --> {dt<dx / Vp}")

Nsteps = int(tmax /dt)

print(f"{λ=} {dx=} {dt=} {Nsteps}")
