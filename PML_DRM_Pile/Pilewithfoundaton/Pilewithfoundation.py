# %%
import pyvista as pv 
import numpy as np
import ctypes 
import os 
import sys
pv.start_xvfb()
# =============================================================================
# information
# =============================================================================
# getting the meshing information from the command line

xwidth, ywidth, zwidth           = (float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]))
eps                              = 1e-6
Xmeshsize, Ymeshsize, Zmeshsize  = (float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]))
PMLThickness                     = np.array([float(sys.argv[7]), float(sys.argv[8]), float(sys.argv[9])])            ; # thickness of the each PML layer
numPMLLayers                     = int(sys.argv[10])                                 ; # number of PML layers
PMLTotalThickness                = PMLThickness * numPMLLayers       ; # total thickness of the PML layers
DRMThickness                     = np.array([float(sys.argv[11]), float(sys.argv[12]), float(sys.argv[13])])            ; # thickness of the DRM layers
numDrmLayers                     = int(sys.argv[14])                                 ; # number of DRM layers
DRMTotalThickness                = DRMThickness * numDrmLayers       ; # total thickness of the DRM layers
padLayers                        = numPMLLayers + numDrmLayers       ; # number of layers to pad the meshgrid
padThickness                     = PMLTotalThickness + DRMThickness  ; # thickness of the padding layers
reg_num_cores                    = int(sys.argv[15])
DRM_num_cores                    = int(sys.argv[16])
PML_num_cores                    = int(sys.argv[17])
Dir                              = sys.argv[18]
OutputDir                        = sys.argv[19]


# xwidth                          = 16
# ywidth                          = 6
# zwidth                          = 8
# eps                             = 1e-6
# Xmeshsize, Ymeshsize, Zmeshsize = (0.5, 0.5, 0.5)
# PMLThickness                    = np.array([Xmeshsize, Ymeshsize, Zmeshsize])            ; # thickness of the each PML layer
# numPMLLayers                    = 2                                 ; # number of PML layers
# PMLTotalThickness               = PMLThickness * numPMLLayers       ; # total thickness of the PML layers
# DRMThickness                    = np.array([Xmeshsize, Ymeshsize, Zmeshsize])            ; # thickness of the DRM layers                                      
# numDrmLayers                    = 1                                 ; # number of DRM layers
# DRMTotalThickness               = DRMThickness * numDrmLayers       ; # total thickness of the DRM layers
# padLayers                       = numPMLLayers + numDrmLayers       ; # number of layers to pad the meshgrid
# padThickness                    = PMLTotalThickness + DRMThickness  ; # thickness of the padding layers
# reg_num_cores                   = 1
# DRM_num_cores                   = 1
# PML_num_cores                   = 3
# Dir                             = "OpenSeesMesh"
# OutputDir                       = "results"

# print the input information
print(f"xwidth: {xwidth}")
print(f"ywidth: {ywidth}")
print(f"zwidth: {zwidth}")
print(f"Xmeshsize: {Xmeshsize}")
print(f"Ymeshsize: {Ymeshsize}")
print(f"Zmeshsize: {Zmeshsize}")
print(f"PMLThickness: {PMLThickness}")
print(f"numPMLLayers: {numPMLLayers}")
print(f"PMLTotalThickness: {PMLTotalThickness}")
print(f"DRMThickness: {DRMThickness}")
print(f"numDrmLayers: {numDrmLayers}")
print(f"DRMTotalThickness: {DRMTotalThickness}")
print(f"padLayers: {padLayers}")
print(f"padThickness: {padThickness}")
print(f"reg_num_cores: {reg_num_cores}")
print(f"DRM_num_cores: {DRM_num_cores}")
print(f"PML_num_cores: {PML_num_cores}")
print(f"Dir: {Dir}")
print(f"OutputDir: {OutputDir}")
# create a dictionary for meshing information
info = {
    "RegularDomain": 1,
    "DRMDomain": 2,
    "PMLDomain": 3,
}
# exit()
# =============================================================================
# define partioner
# =============================================================================
# change the directory to the current directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))

if os.name == 'nt':
    metis_partition_lib = ctypes.CDLL('ll ../../MeshGenerator/lib//Partitioner.dll')
if os.name == 'posix':
    metis_partition_lib = ctypes.CDLL('../../MeshGenerator/lib/libPartitioner.so')

# Define function argument and return types
metis_partition_lib.Partition.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_int32), ctypes.c_int, ctypes.POINTER(ctypes.c_int32)]
metis_partition_lib.Partition.restype = ctypes.c_int

def partition(mesh, numcores):
    numcells = mesh.n_cells
    numpoints = mesh.n_points
    numweights = 1
    cells = np.array(mesh.cells.reshape(-1, 9), dtype=np.int32)
    cells = cells[:,1:]
    cellspointer = cells.flatten().ctypes.data_as(ctypes.POINTER(ctypes.c_int32))
    partitioned = np.empty(numcells, dtype=np.int32)
    partitionedpointer = partitioned.ctypes.data_as(ctypes.POINTER(ctypes.c_int32))
    metis_partition_lib.Partition(numcells, numpoints, cellspointer, numcores, partitionedpointer)
    mesh.cell_data['partitioned'] = partitioned


# =============================================================================
# meshing
# =============================================================================
x = np.arange(-xwidth/2., xwidth/2.+eps, Xmeshsize)
y = np.arange(-ywidth/2., ywidth/2.+eps, Ymeshsize)
z = np.arange(-zwidth, 0+eps, Zmeshsize)



# %% 

# padding x and y for PML and DRM layers
x  = np.pad(x, (numDrmLayers,numDrmLayers), "linear_ramp", end_values=(x.min()-DRMTotalThickness[0], x.max()+DRMTotalThickness[0]))
y  = np.pad(y, (numDrmLayers,numDrmLayers), "linear_ramp", end_values=(y.min()-DRMTotalThickness[1], y.max()+DRMTotalThickness[1]))
z  = np.pad(z, (numDrmLayers,0), "linear_ramp", end_values=(z.min()-DRMTotalThickness[2]))

# padding the x and y for PML and PML layers
x  = np.pad(x, (numPMLLayers,numPMLLayers), "linear_ramp", end_values=(x.min()-PMLTotalThickness[0], x.max()+PMLTotalThickness[0]))
y  = np.pad(y, (numPMLLayers,numPMLLayers), "linear_ramp", end_values=(y.min()-PMLTotalThickness[1], y.max()+PMLTotalThickness[1]))
z  = np.pad(z, (numPMLLayers,0), "linear_ramp", end_values=(z.min()-PMLTotalThickness[2]))

# %%
x, y, z = np.meshgrid(x, y, z, indexing='ij')

mesh = pv.StructuredGrid(x, y, z)

x = np.arange(-1.5, 1.5+eps, Xmeshsize)
y = np.arange(-1.5, 1.5+eps, Ymeshsize)
z = np.arange(0, 1+eps, Zmeshsize/2.0)
mesh.cell_data["matTag"] = np.ones(mesh.n_cells,dtype=np.uint8)

x, y, z = np.meshgrid(x, y, z, indexing='ij')
found = pv.StructuredGrid(x, y, z)
found.cell_data["matTag"] = np.ones(found.n_cells,dtype=np.uint8)*2

mesh =mesh.merge(found,merge_points=True,tolerance=1e-6,progress_bar = True)


# %%
# =============================================================================
# sperate PML layer 
# =============================================================================
bounds = mesh.bounds
xmin = bounds[0] + PMLTotalThickness[0]
xmax = bounds[1] - PMLTotalThickness[0]
ymin = bounds[2] + PMLTotalThickness[1]
ymax = bounds[3] - PMLTotalThickness[1]
zmin = bounds[4] + PMLTotalThickness[2]
zmax = bounds[5] + PMLTotalThickness[2]
cube = pv.Cube(bounds=[xmin,xmax,ymin,ymax,zmin,zmax])
PML = mesh.clip_box(cube,invert=True,crinkle=True,progress_bar = True)
reg = mesh.clip_box(cube,invert=False,crinkle=True,progress_bar = True)
# %%
# now find DRM layer
indices = reg.find_cells_within_bounds([xmin + DRMTotalThickness[0] + eps,
                              xmax - DRMTotalThickness[0] - eps,
                              ymin + DRMTotalThickness[1] + eps,
                              ymax - DRMTotalThickness[1] - eps,
                              zmin + DRMTotalThickness[2] + eps,
                              zmax + DRMTotalThickness[2] + eps])

# now create complemntary indices for DRM
DRMindices = np.ones(reg.n_cells,dtype=bool)
DRMindices[indices] = False
DRMindices = np.where(DRMindices)[0]



reg.cell_data['Domain'] = np.ones(reg.n_cells,dtype=np.int8)*info["DRMDomain"]
reg.cell_data['Domain'][indices] = info["RegularDomain"]
PML.cell_data['Domain'] = np.ones(PML.n_cells,dtype=np.int8)*info["PMLDomain"]
reg.cell_data['partitioned'] = np.zeros(reg.n_cells,dtype=np.int32)


# partitioning regular mesh
regular = reg.extract_cells(indices,progress_bar=True)
DRM     = reg.extract_cells(DRMindices,progress_bar=True)

if reg_num_cores > 1:
    partition(regular,reg_num_cores)
if DRM_num_cores > 1:
    partition(DRM,DRM_num_cores)
if PML_num_cores > 1:
    partition(PML,PML_num_cores)

reg.cell_data['partitioned'][regular["vtkOriginalCellIds"]] = regular.cell_data['partitioned']
reg.cell_data['partitioned'][DRM["vtkOriginalCellIds"]] = DRM.cell_data['partitioned'] + reg_num_cores
PML.cell_data['partitioned'] = PML.cell_data['partitioned'] + reg_num_cores + DRM_num_cores
PML.cell_data["matTag"] = np.ones(PML.n_cells,dtype=np.uint8)
# %%
# merging PML and regular mesh to create a single mesh
mesh = reg.merge(PML,merge_points=False,tolerance=1e-6,progress_bar = True)


# mapping between PML and regular mesh on the boundary
mapping = mesh.clean(produce_merge_map=True)["PointMergeMap"]
regindicies = np.where(mapping[PML.n_points:]<PML.n_points)[0]
PMLindicies = mapping[PML.n_points+regindicies]


mesh.point_data["boundary"] = np.zeros(mesh.n_points,dtype=int)-1
mesh.point_data["boundary"][PMLindicies] = regindicies + PML.n_points
mesh.point_data["boundary"][PML.n_points + regindicies] = PMLindicies 

indices = np.where(mesh.point_data["boundary"]>0)[0]
# %%

#  =============================================================================
# write the mesh
# =============================================================================
if not os.path.exists(Dir):
    os.makedirs(Dir)
else :
    # remove the files in the directory
    for file in os.listdir(Dir):
        os.remove(os.path.join(Dir,file))


min_core = mesh.cell_data['partitioned'].min()
max_core = mesh.cell_data['partitioned'].max()

# write the  mesh nodes
for core  in range(min_core,max_core+1):
    tmp  = mesh.extract_cells(np.where(mesh.cell_data['partitioned']==core)[0])
    f  = open(Dir + "/Nodes" + str(core) + ".tcl", "w")

    for i in range(tmp.n_points):
        f.write(f"node  {tmp['vtkOriginalPointIds'][i]} {tmp.points[i][0]} {tmp.points[i][1]} {tmp.points[i][2]}\n")
    f.close()
# %%
# writing the mesh elements
for core in range(min_core,max_core+1):
    tmp  = mesh.extract_cells(np.where(mesh.cell_data['partitioned']==core)[0])
    f    = open(Dir + "/Elements" + str(core) + ".tcl", "w")
    for eletag in range(tmp.n_cells):
        f.write(f"eval \"element $elementType {tmp['vtkOriginalCellIds'][eletag]} {' '.join(str(x) for x in tmp['vtkOriginalPointIds'][tmp.get_cell(eletag).point_ids])} $matTag{tmp['matTag'][eletag]}\" \n")
    f.close()

# writing the boundary files
for core in range(reg_num_cores + DRM_num_cores , max_core+1):
    tmp = mesh.extract_cells(np.where(mesh.cell_data['partitioned']==core)[0])
    f = open(Dir + "/Boundary" + str(core) + ".tcl", "w")
    for i in range(tmp.n_points):
        if tmp["boundary"][i] != -1:
            x,y,z = tmp.points[i]
            nodeTag1 = tmp['vtkOriginalPointIds'][i]
            nodeTag2 = tmp['boundary'][i]
            f.write(f"node {nodeTag2} {str(x)} {str(y)} {str(z)}\n")
            f.write(f"equalDOF {nodeTag2} {nodeTag1} 1 2 3\n")



# =============================================================================
# printing information
# =============================================================================
print(f"Number of regular cores: {reg_num_cores}")
print(f"Number of DRM cores: {DRM_num_cores}")
print(f"Number of PML cores: {PML_num_cores}")
print(f"Number of regular elements: {regular.n_cells} roughly {int(regular.n_cells/reg_num_cores)} each core")
print(f"Number of DRM elements: {DRM.n_cells} roughly {int(DRM.n_cells/DRM_num_cores)} each core")
print(f"Number of PML elements: {PML.n_cells} roughly {int(PML.n_cells/PML_num_cores)} each core")
print(f"Number of total elements: {mesh.n_cells}")
print(f"Number of total points: {mesh.n_points}")
print(f"Number of cores: {max_core-min_core+1}")

# %%
# mesh.plot(scalars="partitioned",show_edges=True,show_grid=True,show_axes=True,show_bounds=True)
pl = pv.Plotter()
pl.add_mesh(mesh,scalars="matTag",show_edges=True)
pl.export_html(os.path.join(OutputDir,"mesh.html"))
# pl.show()
pl.close()

# 
if not os.path.exists(OutputDir):
    os.makedirs(OutputDir)

mesh.save(os.path.join(OutputDir,"mesh.vtk"),binary=True)




# %%
# chek for each element that the nodes go counter clockwise and first the below surface nodes and then the above surface nodes
# def iscounterclockwise(points):
#     # points = np.array(points)
#     points = points - points.mean(axis=0)
#     signed_area = 0
#     for j in range(points.shape[0]):
#         p1 = points[j,:]
#         p2 = points[(j+1)%points.shape[0],:]
#         x1,y1,_ = p1
#         x2,y2,_ = p2
#         signed_area += (x1*y2 - x2*y1)
    
#     return signed_area/2.0

# for i in range(mesh.n_cells):
#     cell = mesh.get_cell(i)
#     points = mesh.points[cell.point_ids]
    
#      # first surface nodes
#     S1 = points[:4,:]
#     S2 = points[4:,:]
    
#     # check that the S2 is above S1
#     if S1[:,2].mean() > S2[:,2].mean():
#         print(f"Element {i} S1 is above S2")
#     # else :
#     #     print(f"Element {i} S2 is above S1")
        
#     # calcuting signed area of the S1
#     S1 = iscounterclockwise(S1)
#     S2 = iscounterclockwise(S2)
#     if S1 < 0:
#         print(f"Element {i} S1 is not counter clockwise")
#     if S2 < 0:
#         print(f"Element {i} S2 is not counter clockwise")



        


