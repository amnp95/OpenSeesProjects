# %%
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt 
import os
# change the path to the directory where the file is located
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# ==================================================================
# load nodes and elements dataframes
# ==================================================================
mesh = pv.read("results/mesh.vtk")
info = {
    "RegularDomain": 1,
    "DRMDomain": 2,
    "PMLDomain": 3,
}
# filter cells with Domin tag 1
indices = mesh['Domain'] == 1
grid = mesh.extract_cells(indices)
cleangrid = grid.copy()
cleangrid.clear_data()

# ==================================================================
# create interface mesh
# ==================================================================
data = np.loadtxt("results/interfaceInfo.dat")
interface = pv.PolyData(data[:,1:4])

# ==================================================================
# beam mesh
# ==================================================================
beamnodes = np.loadtxt("results/beamNodes.dat")
beamcells = np.loadtxt("results/beamElems.dat", dtype=int)
beamnodes = beamnodes[:,1:]
beamcells[:,0] = 2
beamcells[:,1:] = beamcells[:,1:] - 1
celltypes  = np.ones(beamcells.shape[0],dtype= int) * pv.CellType.LINE
beam = pv.UnstructuredGrid(beamcells.tolist(),celltypes.tolist(),beamnodes.tolist())

# %%
data1 = np.loadtxt("results/NodeAcclPML0.out")
data2 = np.loadtxt("results/NodeAcclPML10.out")
i = grid.find_closest_point([0,0,0])
plt.plot(data1[:,0],(data1[:,1::3])[:,i])
plt.plot(data2[:,0],(data2[:,1::3])[:,i])
plt.legend(["PML=5layers","PML=4layers"])

# %%
data1 = np.loadtxt("results/BeamDispPML0.out")
data2 = np.loadtxt("results/BeamDispPML10.out")

i = beam.find_closest_point([0,0,2])
plt.plot(data1[:,0],(data1[:,1::3])[:,i])
plt.plot(data2[:,0],(data2[:,1::3])[:,i])





# %%
# ==================================================================
# load displacement
# ==================================================================
griddisp   = np.loadtxt("results/NodeDispPML0.out")
beamdisp   = np.loadtxt("results/BeamDispPML0.out")
gridaccel  = np.loadtxt("results/NodeAcclPML0.out")
# griddisp   = np.loadtxt("results/NodeDisp0.out")
# beamdisp   = np.loadtxt("results/BeamDisp0.out")
# gridaccel  = np.loadtxt("results/NodeAccl0.out")
# gridStress = np.loadtxt("results/ElementStress0.out")
# %%
time      = griddisp[:,0]
dt        = time[1] - time[0]
frameRate = 1/dt
griddisp  = griddisp[:,1:]
gridaccel = gridaccel[:,1:]
beamdisp  = beamdisp[:,1:]
griddispx = griddisp[:,0::3]
griddispy = griddisp[:,1::3]
griddispz = griddisp[:,2::3]
dispfactor = 50
print(dt)
print(frameRate)
# %%
slicer = pv.Cube(center=[0,-0.25,0], x_length=150, y_length=0.24, z_length=150)
slicedmesh = cleangrid.copy()
slicedmesh = slicedmesh.clip_box(slicer,invert=False,crinkle=True)
slicepointindexes = slicedmesh["vtkOriginalPointIds"]
slicecelliindexes = slicedmesh["vtkOriginalCellIds"]


direcotry = "./movie/Pics"
if not os.path.exists(direcotry):
    if not os.path.exists("./movie"):
        os.mkdir("./movie")
    os.mkdir(direcotry)
originalcoords     = slicedmesh.points
beamoriginalcoords = beam.points
name  =  "acceleration"
array =  np.abs(gridaccel[:,0::3])

mx = array.max()
mn = array.min()
mn,mx = np.sort(np.abs([mx,mn]))
array = array[:,slicepointindexes]
# empty the directory
files = os.listdir(direcotry)
for file in files:
    os.remove(f"{direcotry}/{file}")

pl = pv.Plotter(off_screen=True)     
mx = 1.0
mn = 0.0
frameRate = 60
# for i in range(500,800,1):
for i in range(0,griddisp.shape[0],1):
    slicedmesh.points = originalcoords     + dispfactor * griddisp[i,:].reshape(-1,3)[slicepointindexes]
    beam.points       = beamoriginalcoords + dispfactor * beamdisp[i,:].reshape(-1,3)
    slicedmesh.point_data[f"{name}"] = array[i,:]
    pl.add_mesh(slicedmesh,scalars=f"{name}", show_edges=True, cmap="coolwarm", clim=[mn,mx])
    pl.add_mesh(beam, color="blue", line_width=5)
    pl.camera_position = 'xz'
    pl.save_graphic(f"{direcotry}/{i:04d}.svg")
    pl.clear()
pl.close()

# %%
direcotry = "./movie/Pics"
frameRate = 60
os.system(f"ffmpeg -y -framerate {frameRate} -pattern_type glob -i '{direcotry}/*.svg' -c:v libx264 -pix_fmt yuv420p tmp.mp4")



# # ==================================================================
# cube = pv.Cube(center=[0,0,0], x_length=150, y_length=0.75, z_length=150)
# gridslice = cleangrid.copy()
# gridslice = gridslice.clip_box(cube,invert=False,crinkle=True)
# sliceindexes = gridslice["vtkOriginalCellIds"]


# pl = pv.Plotter(off_screen=True)
# mx = gridaccel.max()
# mn = gridaccel.min()
# mx = np.abs(mx) if np.abs(mx) > np.abs(mn) else np.abs(mn)
# mn = 0
# dispfactor = 100
# print(gridaccel.shape)

# for i in range(0,griddisp.shape[0],1):
#     grid.points = originalcoords + dispfactor * griddisp[i,:].reshape(-1,3)
#     grid.point_data["acceleration"] = gridaccel[i,:]
#     beam.points = beamoriginalcoords + dispfactor * beamdisp[i,:].reshape(-1,3)
#     pl.add_mesh(grid.clip_box(cube,invert=False,crinkle=True),scalars="acceleration", show_edges=True, cmap="coolwarm", clim=[mn,mx])
#     pl.add_mesh(beam, color="blue", line_width=5)
#     # pl.camera_position = 'xz'
#     pl.save_graphic("./movie/Disp/U/{:04d}.svg".format(i))
#     pl.clear()
# pl.close()
# %%
# %%
# %
# # ==================================================================

# # %%
# # # ==================================================================
# # griddisp   = np.loadtxt("results/NodeDisp_PML.out")
# # time = griddisp[:,0]
# # griddisp = griddisp[:,1:]
# # griddispx = griddisp[:,0::3]
# # index = grid.find_closest_point([0,0,0])
# # index2 = grid.find_closest_point([1,0,0])
# # index3 = grid.find_closest_point([2.0,0,0])
# # plt.plot(time, griddispx[:,index])
# # plt.plot(time, griddispx[:,index2])
# # plt.plot(time, griddispx[:,index3])
# # # beamdispx = beamdisp[:,0::3]
# # # index = beam.find_closest_point([0,0,0])
# # # index2 = beam.find_closest_point([0,0,-1])
# # # index3 = beam.find_closest_point([0,0,-3])
# # # plt.plot(time, beamdispx[:,index])
# # # plt.plot(time, beamdispx[:,index2])
# # # plt.plot(time, beamdispx[:,index3])
# # # plt.xlim(right=1.0)
# # plt.show()


# # # %%
# # gridaccelx = gridaccel[:,0::3]
# # tmpdis = gridaccelx[:,index]


# # # plot the frequncy content of of the tmpdis

# # frequncy = np.fft.fftfreq(tmpdis.size, d=0.001)
# # fft = np.fft.fft(tmpdis)
# # plt.plot(frequncy, np.abs(fft))
# # plt.xlim(0,20)
# # frequncy
# # # %%
# # # Generate a sample signal
# # fs = 1000  # Sampling frequency in Hz
# # t = np.arange(0, 1, 1/fs)  # Time vector of 1 second

# # # Create a signal with two frequencies: 50 Hz and 120 Hz
# # f1 = 50  # Frequency of the first sine wave in Hz
# # f2 = 120  # Frequency of the second sine wave in Hz
# # signal = 0.6 * np.sin(2*np.pi*f1*t) + 0.4 * np.sin(2*np.pi*f2*t)

# # # Perform the FFT
# # fft_result = np.fft.fft(signal)
# # fft_freq = np.fft.fftfreq(len(signal), 1/fs)

# # # Plot the signal
# # plt.figure(figsize=(12, 6))
# # plt.subplot(2, 1, 1)
# # plt.plot(t, signal)
# # plt.title('Time Domain Signal')
# # plt.xlabel('Time (seconds)')
# # plt.ylabel('Amplitude')

# # # Plot the magnitude of the FFT result
# # plt.subplot(2, 1, 2)
# # plt.stem(fft_freq, np.abs(fft_result), 'b', markerfmt=" ", basefmt="-b")
# # plt.title('Frequency Domain Signal')
# # plt.xlabel('Frequency (Hz)')
# # plt.ylabel('Magnitude')
# # plt.xlim(0, fs/2)  # Only plot the positive frequencies

# # plt.tight_layout()
# # plt.show()


# # # %%
# # # ==================================================================
# # # create a movie
# # # ==================================================================
# # plotter = pv.Plotter(off_screen=True)
# # plotter.open_movie("movie.mp4")
# # plotter.add_mesh(grid)
# # originalcoords = grid.points
# # plotter.show(auto_close=False)  # only necessary for an off-screen movie
# # plotter.write_frame()  # write initial data

# # # Update scalars on each frame
# # for i in range(0,1000,5):
# #     grid.points = grid.points + 100 * griddisp[i,:].reshape(-1,3)
# #     plotter.add_text(f"Iteration: {i}", name='time-label')
# #     plotter.write_frame()  # Write this frame

# # # Be sure to close the plotter when finished
# # plotter.close()





# # # %%

# %%
