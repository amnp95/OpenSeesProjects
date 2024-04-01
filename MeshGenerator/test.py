# we want to import D:\Projects\Mesh\Partitioner\build\Release\PartitionLib.dll which is  C library in python and use it 
 # %%
import ctypes
import numpy as np
import os
import pyvista as pv

#// load the library
# find that the system is ubuntu or windows
if os.name == 'nt':
    metis_partition_lib = ctypes.CDLL('./lib/Partitioner.dll')
if os.name == 'posix':
    metis_partition_lib = ctypes.CDLL('./lib/libPartitioner.so')



# Define function argument and return types
metis_partition_lib.Partition.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_int32), ctypes.c_int, ctypes.POINTER(ctypes.c_int32)]
metis_partition_lib.Partition.restype = ctypes.c_int

# %%
# Create a simple mesh
x = np.arange(-20, 20+0.01, 1.0)
y = np.arange(-20, 20+0.01, 1.0)
z = np.arange(-20, 20+0.01, 1.0)

x, y, z = np.meshgrid(x, y, z)

grid = pv.StructuredGrid(x, y, z)

# change grid to unstructured grid
grid = grid.cast_to_unstructured_grid()


numcores = 4
numcells = grid.n_cells
numpoints = grid.n_points
numweights = 1


cells = np.array(grid.cells.reshape(-1, 9), dtype=np.int32)
cells = cells[:,1:]
# %%
cellspointer = cells.flatten().ctypes.data_as(ctypes.POINTER(ctypes.c_int32))
partitioned = np.empty(numcells, dtype=np.int32)
partitionedpointer = partitioned.ctypes.data_as(ctypes.POINTER(ctypes.c_int32))
metis_partition_lib.Partition(numcells, numpoints, cellspointer, numcores, partitionedpointer)
grid.cell_data['partitioned'] = partitioned
grid.plot(show_edges=False, scalars='partitioned', cmap='rainbow')

# %%
