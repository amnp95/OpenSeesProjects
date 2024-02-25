import os
import ctypes
import numpy as np

if os.name == 'nt':
    metis_partition_lib = ctypes.CDLL('./lib/Partitioner.dll')
if os.name == 'posix':
    metis_partition_lib = ctypes.CDLL('./lib/libPartitioner.so')

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

