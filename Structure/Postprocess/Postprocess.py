# %%
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import os

typename = "FIXED"
# ==================================================================
# Structure mesh
# ==================================================================
beamnodes = np.loadtxt(f"results/{typename}/StructureNodes.dat")
beamcells = np.loadtxt(f"results/{typename}/StructureElems.dat", dtype=int)
nodetags  = beamnodes[:,0]
nodetags = nodetags.astype(int)
# create mapping between nodetags the index of the tag in the nodetags arrays
nodetags = dict(zip(nodetags,range(nodetags.shape[0])))

beamnodes = beamnodes[:,1:]


beamcells[:,0]  = 2
beamcells[:,1:] = np.vectorize(nodetags.get)(beamcells[:,1:])
celltypes  = np.ones(beamcells.shape[0],dtype= int) * pv.CellType.LINE
beam = pv.UnstructuredGrid(beamcells.tolist(),celltypes.tolist(),beamnodes.tolist())
# ==================================================================
# load nodes and elements dataframes
# ==================================================================
mesh = pv.read("results/mesh.vtk")
info = {
    "RegularDomain": 1,
    "DRMDomain": 2,
    "PMLDomain": 3,
}
# typename = "FIXED"
typename = "PML"
# filter cells with Domin tag 1
indices = mesh['Domain'] == info["RegularDomain"]
grid = mesh.extract_cells(indices)
cleangrid = grid.copy()
cleangrid.clear_data()


pl = pv.Plotter()
pl.add_mesh(grid, color="white",opacity=0.5,show_edges=True)
pl.add_mesh(beam, color="blue", line_width=3)
pl.show()
# %%
