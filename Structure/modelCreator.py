# %%
import numpy as np
from MeshGenerator.MeshGenrator import *
from MeshGenerator.Infowriter import *
import os
import shutil

# ============================================================================
# Cores information
# ============================================================================
regcores       = 5
pmlcores       = 3
drmcores       = 1
structurecores = 1
AnalysisType   = "PML"                            # options: "DRMFIXED", PMLDRM", "ASDADRM"
Target         = "Soil-with-structure"            # options: "Soil", "Structure", "Soil-with-structure"

# ============================================================================
# Structure information
# ============================================================================
StructureType = "STEEL"                             # Options: STEEL, CONCRETE, Custom
NStory = 10                                         # number of stories above ground level   
NBay   = 4                                          # number of bays in X direction
NBayZ  = 4                                          # number of bays in Y direction
StartNodeX = -9.0                                   # X coordinate of the first node
StartNodeY = -9.0                                   # Y coordinate of the first node
StartNodeZ = -1.5                                    # Z coordinate of the first node
meter = 1.0                                         # meter to specified unit conversion (SI unit)
ft    = 0.3048                                      # feet to meter conversion (SI unit)
LCol  = 3*meter                                     # column height (parallel to Z axis)
LBeam = 4.5*meter                                   # beam length (parallel to X axis)
LGird = 4.5*meter                                   # girder length (parallel to Y axis)
SectionType = "Elastic"                             # options: Elastic, FiberSection
HaveStructure = "YES"                               # options: "YES", "NO"

# ============================================================================
# Soil information
# ============================================================================
dx            = 2.5
dy            = dx
dz            = dx
llx           = 60.0
lly           = 60.0
llz           = 40.0
drmthicknessx = dx
drmthicknessy = dy
drmthicknessz = dz
numdrmlayers  = 2
lx            = llx - 2*numdrmlayers*drmthicknessx
ly            = lly - 2*numdrmlayers*drmthicknessy
lz            = llz - 1*numdrmlayers*drmthicknessz
nx            = lx/dx
ny            = ly/dy
nz            = lz/dz


# ============================================================================
# PML information
# ============================================================================
AbsorbingElements = "PML"                           # could be "ASDA" or "NORMAL"
numpmllayers  = 2
pmlthicknessx = dx
pmlthicknessy = dy
pmlthicknessz = dz
pmltotalthickness = numpmllayers*pmlthicknessx
HaveAbsorbingElements = "NO"                      # options: "YES", "NO"


# ============================================================================
# General information
# ============================================================================
meshdir       = f"OpenSeesmesh/{AnalysisType}"
outputdir     = f"Results/{AnalysisType}"
DRMFile       = "/home/amnp95/Projects/OpenSeesProjects/Structure/DRMload/SurfaceWave.h5drm"

# ============================================================================
# Embedding foundation
# ============================================================================
EmbeddingFoundation = "YES" # options: "YES", "NO"
EmbeddedFoundation = {
    "xmax":  10.0,
    "xmin": -10.0,
    "ymax":  10.0,
    "ymin": -10.0,
    "zmax":  0.0,
    "zmin": -5.0,
}

# ============================================================================
# Fondation information
# ============================================================================
HaveFoundation = "YES"
AttachFoundation = "NO"
foundationBlocks = []

# adding a block
foundationBlocks.append({
    "matTag": 1,
    "xmax":  10.0,
    "xmin": -10.0,
    "ymax":  10.0,
    "ymin": -10.0,
    "zmax":  -1.5,
    "zmin":  -4.5,
    "Xmeshsize" : 1.0,
    "Ymeshsize" : 1.0,
    "Zmeshsize" : 1.0,
})


# ============================================================================
# piles information
# ============================================================================
pilelist = []

x = np.arange(-7.0, 7.0+1e-6, 3.5)
y = np.arange(-7.0, 7.0+1e-6, 3.5)

x, y = np.meshgrid(x, y)
x = x.flatten()
y = y.flatten()

HavePiles = "YES"
for i in range(len(x)):
    pilelist.append({
        "xtop":  x[i],
        "ytop":  y[i],
        "ztop": -3.0,
        "xbottom":  x[i],
        "ybottom":  y[i],
        "zbottom": -10.0,
        "numberofElements": 6,
    })

# ============================================================================
# Create directories
# ============================================================================
# delete the directories if exists
if os.path.exists(meshdir):
    # removing directory
    shutil.rmtree(meshdir, ignore_errors=False)

if os.path.exists(outputdir):
    shutil.rmtree(outputdir, ignore_errors=False) 


# create directories first 
if not os.path.exists(meshdir):
    dirs = meshdir.split("/")
    for i in range(1, len(dirs)+1):
        if not os.path.exists("/".join(dirs[:i])):
            os.makedirs("/".join(dirs[:i]))

if not os.path.exists(outputdir):
    dirs = outputdir.split("/")
    for i in range(1, len(dirs)+1):
        if not os.path.exists("/".join(dirs[:i])):
            os.makedirs("/".join(dirs[:i]))



# ============================================================================
# Creating the mesh
# ============================================================================
info = {
    "xwidth": lx,
    "ywidth": ly,
    "zwidth": lz,
    "Xmeshsize": dx,
    "Ymeshsize": dy,
    "Zmeshsize": dz,
    "PMLThickness": [pmlthicknessx, pmlthicknessy, pmlthicknessz],
    "numPMLLayers": numpmllayers,
    "DRMThickness": [drmthicknessx, drmthicknessy, drmthicknessz],
    "numDrmLayers": numdrmlayers,
    "reg_num_cores": regcores,
    "DRM_num_cores": drmcores,
    "PML_num_cores": pmlcores,
    "Structure_num_cores": structurecores,
    "Dir": meshdir,
    "OutputDir": outputdir,
    "AbsorbingElements": AbsorbingElements,
    "DRMfile": DRMFile,
    "EmbeddingFoundation": EmbeddingFoundation,
    "EmbeddedFoundationDict": EmbeddedFoundation,
    "HaveFoundation": HaveFoundation,
    "foundationBlocks": foundationBlocks,
    "pilelist": pilelist,
    "HavePiles": HavePiles,
    "HaveStructure": HaveStructure,
    "AttachFoundation": AttachFoundation,

}

numcells, numpoints = DRM_PML_Foundation_Meshgenrator(info)
print(f"Number of cells: {numcells}")
print(f"Number of points: {numpoints}")
# ============================================================================
# Writing the information file
# ============================================================================
info = {
    "soilfoundation_num_cells": numcells,
    "soilfoundation_num_points": numpoints,
    "AnalysisType": AnalysisType,
    "regcores": regcores,
    "pmlcores": pmlcores,
    "drmcores": drmcores,
    "structurecores": structurecores,
    "StructureType": StructureType,
    "NStory": NStory,
    "NBay": NBay,
    "NBayZ": NBayZ,
    "StartNodeX": StartNodeX,
    "StartNodeY": StartNodeY,
    "StartNodeZ": StartNodeZ,
    "meter": meter,
    "ft": ft,
    "LCol": LCol,
    "LBeam": LBeam,
    "LGird": LGird,
    "SectionType": SectionType,
    "dx": dx,
    "dy": dy,
    "dz": dz,
    "llx": llx,
    "lly": lly,
    "llz": llz,
    "drmthicknessx": drmthicknessx,
    "drmthicknessy": drmthicknessy,
    "drmthicknessz": drmthicknessz,
    "numdrmlayers": numdrmlayers,
    "lx": lx,
    "ly": ly,
    "lz": lz,
    "nx": nx,
    "ny": ny,
    "nz": nz,
    "AbsorbingElements": AbsorbingElements,
    "numpmllayers": numpmllayers,
    "pmlthicknessx": pmlthicknessx,
    "pmlthicknessy": pmlthicknessy,
    "pmlthicknessz": pmlthicknessz,
    "pmltotalthickness": pmltotalthickness,
    "HaveAbsorbingElements": HaveAbsorbingElements,
    "meshdir": meshdir,
    "outputdir": outputdir,
    "DRMFile": DRMFile,
    "EmbeddingFoundation": EmbeddingFoundation,
    "EmbeddedFoundation": EmbeddedFoundation,
    "HaveFoundation": HaveFoundation,
    "foundationBlocks": foundationBlocks,
    "pilelist": pilelist,
    "HavePiles": HavePiles,
    "HaveStructure": HaveStructure,
}
infowriter(info,meshdir)

# ============================================================================
# copy the related file as model.tcl t the current directory
# ============================================================================
def copy_file(source_path, destination_path):
    with open(destination_path, 'wb') as dst_file:
        with open(f"{meshdir}/Modelinfo.tcl", 'rb') as src_file:
            dst_file.write(src_file.read())
        with open(source_path, 'rb') as src_file:
            dst_file.write(src_file.read())


# delete the model file if exists
if os.path.exists("./model.tcl"):
    os.remove("./model.tcl")

if Target == "Soil-with-structure":
    copy_file(f"MeshGenerator/models/Soil_with_structure.tcl", f"./model.tcl")

# %%
