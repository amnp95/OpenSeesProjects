# %%
import numpy as np
import pyvista as pv
import json
import os
from Postprocess.Postproceesfunctions import *
ResultsPath = "Results/PML"

# change the path to the directory where the file is located
path = os.path.dirname(os.path.abspath(__file__))
os.chdir(path)
# ==================================================================
# load the modelInfo.json file
# ==================================================================
json_file = open(f'{ResultsPath}/ModelInfo.json')
model_info = json.load(json_file)

json_file.close()
info = {
    "Foundation": 1,
    "RegularDomain": 2,
    "DRMDomain": 3,
    "PMLDomain": 4,
}
# %%
# ==================================================================
# Structure mesh
# ==================================================================
StructureNodes = np.loadtxt(f"{ResultsPath}/StructureNodes.dat")
StructureElems = np.loadtxt(f"{ResultsPath}/StructureElements.dat", dtype=int)
nodetags  = StructureNodes[:,0]
nodetags = nodetags.astype(int)
# create mapping between nodetags the index of the tag in the nodetags arrays
nodetags = dict(zip(nodetags,range(nodetags.shape[0])))
StructureNodes = StructureNodes[:,1:]
StructureElems[:,0] = 2
StructureElems[:,1:] = np.vectorize(nodetags.get)(StructureElems[:,1:])
celltypes  = np.ones(StructureElems.shape[0],dtype= int) * pv.CellType.LINE
Structure = pv.UnstructuredGrid(StructureElems.tolist(),celltypes.tolist(),StructureNodes.tolist())
# ==================================================================
# structure foundation interface beams
# ==================================================================
StructureFoundationNodes = np.loadtxt(f"{ResultsPath}/StructureFoundationBeamNodes.dat")
StructureFoundationElems = np.loadtxt(f"{ResultsPath}/StructureFoundationBeamElements.dat", dtype=int)
nodetags  = StructureFoundationNodes[:,0]
nodetags = nodetags.astype(int)
# create mapping between nodetags the index of the tag in the nodetags arrays
nodetags = dict(zip(nodetags,range(nodetags.shape[0])))
StructureFoundationNodes = StructureFoundationNodes[:,1:]
StructureFoundationElems[:,0] = 2
StructureFoundationElems[:,1:] = np.vectorize(nodetags.get)(StructureFoundationElems[:,1:])
celltypes  = np.ones(StructureFoundationElems.shape[0],dtype= int) * pv.CellType.LINE
StructureFoundation = pv.UnstructuredGrid(StructureFoundationElems.tolist(),celltypes.tolist(),StructureFoundationNodes.tolist())


# ==================================================================
# pile mesh
# ==================================================================
PileNodes = np.loadtxt(f"{ResultsPath}/PileNodes.dat")
PileElems = np.loadtxt(f"{ResultsPath}/PileElements.dat", dtype=int)
nodetags  = PileNodes[:,0]
nodetags = nodetags.astype(int)
# create mapping between nodetags the index of the tag in the nodetags arrays
nodetags = dict(zip(nodetags,range(nodetags.shape[0])))
PileNodes = PileNodes[:,1:]
PileElems[:,0] = 2
PileElems[:,1:] = np.vectorize(nodetags.get)(PileElems[:,1:])
celltypes  = np.ones(PileElems.shape[0],dtype= int) * pv.CellType.LINE
Pile = pv.UnstructuredGrid(PileElems.tolist(),celltypes.tolist(),PileNodes.tolist())
# %%
# soil mesh
mesh  = pv.read(f"{ResultsPath}/mesh.vtk")
indices =  (mesh['Domain'] == info["Foundation"])
pilblock = mesh.extract_cells(indices)
Soil     = mesh.extract_cells((mesh['Domain'] == info["RegularDomain"])) 
pl = pv.Plotter()
Structure["partitioned"] = np.zeros(Structure.n_points)
pl.add_mesh(Structure, color="blue", line_width=4.0,scalars = "partitioned")
# pl.add_mesh(StructureFoundation, color="purple", line_width=4.0)
Pile["partitioned"] = np.zeros(Pile.n_points) +1
Soil["partitioned"] += 1
pilblock["partitioned"] += 1
pilblock["partitioned"] = np.zeros(pilblock.n_points) +1
pl.add_mesh(Pile, color="red", line_width=4.0,scalars = "partitioned")
pl.add_mesh(Soil, color="white",opacity=1.,show_edges=False,scalars = "partitioned")
pl.add_mesh(pilblock, color="grey",opacity=1.,show_edges=False,scalars = "partitioned")

# ==================================================================
# interface points mesh
# ==================================================================
files = os.listdir(ResultsPath)
files = [f for f in files if "PileInterfaceElements" in f]
numpiles = len(files)
for i in range(numpiles):
    InterfaceNodes = np.loadtxt(f"{ResultsPath}/PileInterfaceElements{i}.dat", delimiter="\t",dtype=float)
    InterfaceNodes = InterfaceNodes[:,1:]
    tmpmesh = pv.PolyData(InterfaceNodes)
    tmpmesh["partitioned"] = np.zeros(tmpmesh.n_points)+1
    if i == 0:
        Interface = pv.PolyData(InterfaceNodes)
    if i > 0:
        Interface = Interface.merge(tmpmesh)
    pl.add_mesh(tmpmesh, color="green", point_size=5, render_points_as_spheres=True,interpolate_before_map=True,scalars = "partitioned")

# ==================================================================
# structure foundation interface points
# ==================================================================
files = os.listdir(ResultsPath)
files = [f for f in files if "StructureFoundationInterfaceElements" in f]
numinterfaces = len(files)
for i in range(numinterfaces):
    InterfaceNodes = np.loadtxt(f"{ResultsPath}/StructureFoundationInterfaceElements{i}.dat", delimiter="\t",dtype=float)
    InterfaceNodes = InterfaceNodes[:,1:]
    tmpmesh = pv.PolyData(InterfaceNodes)
    xmin,xmax,ymin,ymax,zmin,zmax = tmpmesh.bounds
    xmean,ymean,zmean= tmpmesh.points.mean(axis=0)
    # create line element with xmean, ymean, zmin and zmax
    line = pv.Line([xmean,ymean,zmin],[xmean,ymean,zmax])
    tmpmesh["partitioned"] = np.zeros(tmpmesh.n_points)+1
    line["partitioned"] = np.zeros(1)+1
    pl.add_mesh(tmpmesh, color="purple", point_size=5, render_points_as_spheres=True,scalars = "partitioned")
    pl.add_mesh(line, color="black", line_width=4.0,scalars = "partitioned")
    


pl.show()


# %%
pileblock = mesh.extract_cells((mesh['partitioned'] == 0))
block = mesh.extract_cells((mesh['Domain'] == 1))
reg = mesh.extract_cells((mesh['Domain'] == 2))
drm = mesh.extract_cells((mesh['Domain'] == 3))
pml = mesh.extract_cells((mesh['Domain'] == 4))
pl = pv.Plotter()
pl.add_mesh(Structure, color="black", line_width=4.0, label="Structure")
# pl.add_mesh(mesh, color="grey",opacity=0.5,show_edges=False,scalars = "Domain", categories=True)
pl.add_mesh(reg, color="royalblue",opacity=.5,show_edges=True, label="Soil Layer")
pl.add_mesh(block, color="grey",opacity=.5,show_edges=False,label="Foundation")
# pl.add_mesh(drm, color="green",opacity=.5,show_edges=False,label="Domain Reduction Layer")
# pl.add_mesh(pml, color="orange",opacity=.5,show_edges=False, label="Perfectly Matched Layer")
pl.add_mesh(Pile, color="red", line_width=4.0)
# ==================================================================
# interface points mesh
# ==================================================================
files = os.listdir(ResultsPath)
files = [f for f in files if "PileInterfaceElements" in f]
numpiles = len(files)
for i in range(numpiles):
    InterfaceNodes = np.loadtxt(f"{ResultsPath}/PileInterfaceElements{i}.dat", delimiter="\t",dtype=float)
    InterfaceNodes = InterfaceNodes[:,1:]
    if i == 0:
        Interface = pv.PolyData(InterfaceNodes)
    if i > 0:
        Interface = Interface.merge(pv.PolyData(InterfaceNodes))

pl.add_mesh(Interface, color="green", point_size=5, render_points_as_spheres=True)

# ==================================================================
# structure foundation interface points
# ==================================================================
files = os.listdir(ResultsPath)
files = [f for f in files if "StructureFoundationInterfaceElements" in f]
numinterfaces = len(files)
for i in range(numinterfaces):
    InterfaceNodes = np.loadtxt(f"{ResultsPath}/StructureFoundationInterfaceElements{i}.dat", delimiter="\t",dtype=float)
    InterfaceNodes = InterfaceNodes[:,1:]
    tmpmesh = pv.PolyData(InterfaceNodes)
    xmin,xmax,ymin,ymax,zmin,zmax = tmpmesh.bounds
    xmean,ymean,zmean= tmpmesh.points.mean(axis=0)
    # create line element with xmean, ymean, zmin and zmax
    line = pv.Line([xmean,ymean,zmin],[xmean,ymean,zmax])
    pl.add_mesh(tmpmesh, color="purple", point_size=5, render_points_as_spheres=True)
    pl.add_mesh(line, color="black", line_width=4.0)
pl.add_legend(bcolor="white",loc="upper left", border=True,size=[0.3,0.3],face="rectangle")
pl.show()
# %%




ResultsPath = "Results/PML"

# ==================================================================
# load the modelInfo.json file
# ==================================================================
json_file = open(f'{ResultsPath}/ModelInfo.json')
model_info = json.load(json_file)

json_file.close()
info = {
    "Foundation": 1,
    "RegularDomain": 2,
    "DRMDomain": 3,
    "PMLDomain": 4,
}

# ==================================================================
# load the Soil mesh 
# ==================================================================
mesh  = pv.read(f"{ResultsPath}/mesh.vtk")
indices = (mesh['Domain'] == info["RegularDomain"]) | (mesh['Domain'] == info["Foundation"])
Soil = mesh.extract_cells(indices)
Gravdisp  = loadPointData(mesh,"GravityNodeDisp", ResultsPath,StartingCore=1,returnSubset=False)
mesh["GravityNodeDisp"] = Gravdisp[-1,:].reshape(-1,3)
mesh.points = mesh.points + Gravdisp[-1,:].reshape(-1,3)*100
pl = pv.Plotter()
pl.add_mesh(mesh, color="white",opacity=1.0,show_edges=True,scalars = "GravityNodeDisp",cmap="coolwarm")
pl.show()


# %%
Gravdisp = loadPointData(mesh,"GravityNodeDisp", ResultsPath,StartingCore=1,returnSubset=True)
Soil["GravityNodeDisp"] = Gravdisp[-1,:].reshape(-1,3)
Soil.points = Soil.points + Gravdisp[-1,:].reshape(-1,3)
pl = pv.Plotter()
pl.add_mesh(Soil, color="white",opacity=1.0,show_edges=True,scalars = "GravityNodeDisp",cmap="coolwarm")
pl.add_mesh(Pile, color="red", line_width=4.0)
pl.add_mesh(Structure, color="blue", line_width=4.0)
pl.add_mesh(StructureFoundation, color="purple", line_width=4.0)
pl.show()


# %%
GraveStress  = loadCellData(mesh,"GravityElementStress", ResultsPath, "stdbrick",StartingCore=1,returnSubset=True)
Soil["GravityStress"] = GraveStress[-1,2::6]
bounds = [-100,100,0,100,-100,100]
a
pl = pv.Plotter()
pl.add_mesh(Soil.clip_box(bounds), color="white",opacity=1.0,show_edges=False,scalars = "GravityStress",cmap="coolwarm",interpolate_before_map=True,clim=[a,b])
pl.add_mesh(Pile, color="red", line_width=4.0)
pl.add_mesh(Structure, color="blue", line_width=4.0)
pl.add_mesh(StructureFoundation, color="purple", line_width=4.0)
pl.show()

# %%
# ==================================================================
# load the Soil mesh 
# ==================================================================
mesh  = pv.read(f"{ResultsPath}/mesh.vtk")
indices = (mesh['Domain'] == info["RegularDomain"]) | (mesh['Domain'] == info["Foundation"])
Soil = mesh.extract_cells(indices)
data = np.loadtxt(f"{ResultsPath}/GravityElementStress0.out")
Sx = data[:,0::6].reshape((2,-1,8)).mean(axis=2)
Soil["GravityStress"] = Sx[-1,:]
bounds = [-100,100,0,100,-100,100]
pl = pv.Plotter()
pl.add_mesh(Soil, color="white",opacity="linear",show_edges=False,scalars = "GravityStress",cmap="coolwarm",interpolate_before_map=True)
pl.show()



# %%

# Gravdisp   = loadPointData(mesh,"GravityNodeDisp", ResultsPath,StartingCore=0)
# mesh["GravityNodeDisp"] = Gravdisp[-1,:].reshape(-1,3)
# mesh.points = mesh.points + Gravdisp[-1,:].reshape(-1,3)*100000


# pl = pv.Plotter()
# pl.add_mesh(mesh, color="white",opacity=1.0,show_edges=True,scalars = "GravityNodeDisp",cmap="coolwarm")
# pl.show()
# # ==================================================================
# Gravdisp  = loadPointData(mesh,"GravityNodeDisp", ResultsPath,StartingCore=0,returnSubset=True)
# Soil["GravityNodeDisp"] = Gravdisp[-1,:].reshape(-1,3)
# Soil.points = Soil.points + Gravdisp[-1,:].reshape(-1,3)*100000
# pl = pv.Plotter()
# pl.add_mesh(Soil, color="white",opacity=1.0,show_edges=True,scalars = "GravityNodeDisp",cmap="coolwarm")
# pl.show()

# # ==================================================================
# # load the pile mesh
# # ==================================================================
# pl = pv.Plotter()
# Soil2 = mesh.extract_cells(indices)
# pl.add_mesh(Soil2, color="white",opacity=1.0,show_edges=True)
# pl.show()


#%%




# %%

# ==================================================================
# Structure mesh
# ==================================================================
StructureNodes = np.loadtxt(f"{ResultsPath}/StructureNodes.dat")
StructureElems = np.loadtxt(f"{ResultsPath}/StructureElements.dat", dtype=int)
nodetags  = StructureNodes[:,0]
nodetags = nodetags.astype(int)
# create mapping between nodetags the index of the tag in the nodetags arrays
nodetags = dict(zip(nodetags,range(nodetags.shape[0])))
StructureNodes = StructureNodes[:,1:]
StructureElems[:,0] = 2
StructureElems[:,1:] = np.vectorize(nodetags.get)(StructureElems[:,1:])
celltypes  = np.ones(StructureElems.shape[0],dtype= int) * pv.CellType.LINE
Structure = pv.UnstructuredGrid(StructureElems.tolist(),celltypes.tolist(),StructureNodes.tolist())
# ==================================================================
# structure foundation interface beams
# ==================================================================
StructureFoundationNodes = np.loadtxt(f"{ResultsPath}/StructureFoundationBeamNodes.dat")
StructureFoundationElems = np.loadtxt(f"{ResultsPath}/StructureFoundationBeamElements.dat", dtype=int)
nodetags  = StructureFoundationNodes[:,0]
nodetags = nodetags.astype(int)
# create mapping between nodetags the index of the tag in the nodetags arrays
nodetags = dict(zip(nodetags,range(nodetags.shape[0])))
StructureFoundationNodes = StructureFoundationNodes[:,1:]
StructureFoundationElems[:,0] = 2
StructureFoundationElems[:,1:] = np.vectorize(nodetags.get)(StructureFoundationElems[:,1:])
celltypes  = np.ones(StructureFoundationElems.shape[0],dtype= int) * pv.CellType.LINE
StructureFoundation = pv.UnstructuredGrid(StructureFoundationElems.tolist(),celltypes.tolist(),StructureFoundationNodes.tolist())


# ==================================================================
# pile mesh
# ==================================================================
PileNodes = np.loadtxt(f"{ResultsPath}/PileNodes.dat")
PileElems = np.loadtxt(f"{ResultsPath}/PileElements.dat", dtype=int)
nodetags  = PileNodes[:,0]
nodetags = nodetags.astype(int)
# create mapping between nodetags the index of the tag in the nodetags arrays
nodetags = dict(zip(nodetags,range(nodetags.shape[0])))
PileNodes = PileNodes[:,1:]
PileElems[:,0] = 2
PileElems[:,1:] = np.vectorize(nodetags.get)(PileElems[:,1:])
celltypes  = np.ones(PileElems.shape[0],dtype= int) * pv.CellType.LINE
Pile = pv.UnstructuredGrid(PileElems.tolist(),celltypes.tolist(),PileNodes.tolist())



# %%
pl = pv.Plotter()
pl.add_mesh(Structure, color="blue", line_width=4.0)
pl.add_mesh(StructureFoundation, color="purple", line_width=4.0)
pl.add_mesh(Soil, color="white",opacity=.5,show_edges=False)
pl.add_mesh(Pile, color="red", line_width=4.0)
pl.show()
# %%
import numpy as np
import pyvista as pv
import json
import os
from Postprocess.Postproceesfunctions import *
ResultsPath = "Results/PML"
# ==================================================================
# pile mesh
# ==================================================================
PileNodes = np.loadtxt(f"{ResultsPath}/PileNodes.dat")
PileElems = np.loadtxt(f"{ResultsPath}/PileElements.dat", dtype=int)
nodetags  = PileNodes[:,0]
nodetags = nodetags.astype(int)
# create mapping between nodetags the index of the tag in the nodetags arrays
nodetags = dict(zip(nodetags,range(nodetags.shape[0])))
PileNodes = PileNodes[:,1:]
PileElems[:,0] = 2
PileElems[:,1:] = np.vectorize(nodetags.get)(PileElems[:,1:])
celltypes  = np.ones(PileElems.shape[0],dtype= int) * pv.CellType.LINE
Pile = pv.UnstructuredGrid(PileElems.tolist(),celltypes.tolist(),PileNodes.tolist())

# ==================================================================
# soil mesh
# ==================================================================
Soil  = pv.read(f"{ResultsPath}/mesh.vtk")
# ==================================================================
pl = pv.Plotter()
# pl.add_mesh(Structure, color="blue", line_width=4.0)
pl.add_mesh(Pile, color="red", line_width=4.0)
pl.add_mesh(Soil, color="white",opacity=.5,show_edges=False)
# ==================================================================
# interface points mesh
# ==================================================================
files = os.listdir(ResultsPath)
files = [f for f in files if "PileInterfaceElements" in f]
numpiles = len(files)
for i in range(numpiles):
    InterfaceNodes = np.loadtxt(f"{ResultsPath}/PileInterfaceElements{i}.dat", delimiter="\t",dtype=float)
    InterfaceNodes = InterfaceNodes[:,1:]
    tmpmesh = pv.PolyData(InterfaceNodes)
    if i == 0:
        Interface = pv.PolyData(InterfaceNodes)
    if i > 0:
        Interface = Interface.merge(tmpmesh)
    pl.add_mesh(tmpmesh, color="green", point_size=5, render_points_as_spheres=True,interpolate_before_map=True )
pl.show()
# %%
# ==================================================================
# structure foundation interface points
# ==================================================================
files = os.listdir(ResultsPath)
files = [f for f in files if "StructureFoundationInterfaceElements" in f]
numinterfaces = len(files)
for i in range(numinterfaces):
    InterfaceNodes = np.loadtxt(f"{ResultsPath}/StructureFoundationInterfaceElements{i}.dat", delimiter="\t",dtype=float)
    InterfaceNodes = InterfaceNodes[:,1:]
    tmpmesh = pv.PolyData(InterfaceNodes)
    pl.add_mesh(tmpmesh, color="purple", point_size=5, render_points_as_spheres=True)

# %%
# change the path to the directory where the file is located
import numpy as np
import pyvista as pv
import os
from Postprocess.Postproceesfunctions import *


# ==================================================================
# load nodes and elements dataframes
# ==================================================================
ResultsPath = "Results/PML"
Soil  = pv.read(f"{ResultsPath}/mesh.vtk")
info = {
    "RegularDomain": 1,
    "DRMDomain": 2,
    "PMLDomain": 3,
}
# typename = "FIXED"
typename = "PML"
# filter cells with Domin tag 1
indices = Soil['Domain'] == info["RegularDomain"]
grid = Soil.extract_cells(indices)
cleangrid = grid.copy()
cleangrid.clear_data()

# ==================================================================
# load displacement
# ==================================================================
# typename = "PML"
# typename = "FIXEDNoDamping"
Dir        = f"results/{typename}"
cores      = grid["partitioned"].max()
griddisp   = loadPointData(Soil,"GravityNodeDisp", ResultsPath , cores=cores, timesteps=2, Nodetagsfile="nodeOuputTags",subsetids=grid["vtkOriginalPointIds"])
griddisp.shape[-1]/3.
grid.n_points
# %%
grid["GravityNodeDisp"] = griddisp[1,:].reshape(-1,3)
pl = pv.Plotter()
pl.add_mesh(grid, color="white",opacity=1.0,show_edges=True,scalars = "GravityNodeDisp",cmap="coolwarm")
pl.show()
# %%
