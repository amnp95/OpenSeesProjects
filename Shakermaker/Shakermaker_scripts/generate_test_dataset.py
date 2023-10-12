from shakermaker import shakermaker
from shakermaker.crustmodel import CrustModel
from shakermaker.pointsource import PointSource
from shakermaker.faultsource import FaultSource
from shakermaker.station import Station
from shakermaker.stationlist import StationList
from shakermaker.tools.plotting import ZENTPlot
# from shakermaker.stf_extensions import Discrete
from shakermaker.stf_extensions import SRF2

from math import sqrt
import numpy as np

#Simulation parameters
do_DRM = False
dt = 0.01      # Timestep for output dataset
nfft = 512*2   # Number of samples (need power of 2)
dk = .02       # (wavelength discretization) 
tb = 0         # Zero padding for output (if source-station are too close, needed to avoid aliasing)
tmin = 0.      # start time for simulation
tmax = 60.     # end time for simulation (must be contained in nfft window)

#Single layer crust
crust = CrustModel(2)
Qp=1000.                       # Q-factor for P-wave
Qs=1000.                       # Q-factor for S-wave
Vs = 0.200                     # S-wave speed (km/s)
nu = 0.25                      # Poisson
rho=2.000                      # Mass density (gr/cm**3)
G = rho*Vs**2
M = 2*G*(1-nu)/(1-2*nu)
Vp = sqrt(M/rho)
crust.add_layer(1, Vp, Vs, rho, Qp, Qs)
crust.add_layer(0., Vp, Vs, rho, Qp, Qs)


# Initialize a source time function with
# unit slip, Tr = 1 (rise time) and Tp = 0.1 (peak time)
Tp = 0.1  #peak time
Tr = 1.0  #rise time
stf = SRF2(Tp=Tp, Tr=Tr, dt=dt/10)

#Create a east dip-slip source with 30° dip at 1km depth 
#below the origin
s,d,r = 0., 30., -90.     # Fault plane angles (deg)

for z in [1.]:#[0.2,0.5,1,1.5]:
    source = PointSource([0,0,z], [s,d,r], tt=0, stf=stf)
    fault = FaultSource([source], metadata={"name":"source"})

#Create recording station 30 km away to the east
x0,y0 = 0.,30.           # Station location

#Use distances to compute approximate arrival times for 
#P and S waves (consider a delayes source peak time)
d = np.sqrt(z**2 + x0**2 + y0**2)

tp = d/Vp + Tp
ts = d/Vs + Tp

if not do_DRM:
        s = Station([x0,y0,0], 
                metadata={
                "name":"Your House"
                })
        stations = StationList([s], metadata=s.metadata)

        model = shakermaker.ShakerMaker(crust, fault, stations)
        model.run(
         dt=dt,
         nfft=nfft,
         dk=dk,
         tb=tb,
         tmin=tmin,
         tmax=tmax,
         )

        print(f"{tp=}")
        print(f"{ts=}")

        # Visualize results
        import matplotlib.pyplot as plt
        # ZENTPlot(s, show=True, differentiate=True)#, xlim=[0,3])
        ZENTPlot(s, show=False, integrate=False)#, xlim=[0,3])
        axs = plt.gcf().axes
        for ax in axs:
            ax.axvline(tp,color="r")
            ax.text(tp,0,"P",color="r")
            ax.axvline(ts,color="b")
            ax.text(ts,0,"S",color="b")
        plt.show()
        # ZENTPlot(s, show=True, integrate=True)#, xlim=[0,3])

else:
        from shakermaker.slw_extensions import DRMHDF5StationListWriter
        from shakermaker.sl_extensions import DRMBox
        
        #DRM Box Specification
        _m = 1e-3
        fmax = 20.
        λ = Vs/fmax
        dx = λ / 10
        print(f"{dx=}")
        Lx = 20*_m
        Ly = 20*_m
        Lz = 20*_m
        nx, ny, nz = int(Lx/dx), int(Ly/dx), int(Lz/dx)

        # dx = 10
        # nx, ny, nz = 4, 4, 2
        print(f"{nx=} {ny=} {nz=}")
        x0 = [x0, y0, 0]
        print("Begin DRM define")
        stations = DRMBox(x0,
            [nx,ny,nz],
            [dx,dx,dx],
            metadata={"name":"test"})

 
        #H5DRM writer
        writer = DRMHDF5StationListWriter("test.h5drm")



        #Create model, set parameters and run
        model = shakermaker.ShakerMaker(crust, fault, stations)
        model.run(
         dt=dt,
         nfft=nfft,
         dk=dk,
         tb=tb,
         tmin=tmin,
         tmax=tmax,
         writer=writer
         )
        # 20.48

print("Done!")
exit(0)
