import pyvista as pv
import numpy as np


def Domain(nodes, elements):
    '''
    View soil and PML domains of the 2D box model.

    inputs
    ------------
    nodes : pandas.DataFrame 
        The nodes of the mesh.
    elements : pandas.DataFrame
        The elements of the mesh.
    
        
    outputs
    ------------
    None
    '''
    # initialize plotter and its color cycler
    pl = pv.Plotter()
    pl.set_color_cycler("default")

    # filter elements
    for dm in ["Soil", "PML", "Structure"]:
        soileles = elements[elements['Domain'] == dm]
        soilnodes = nodes[nodes['Domain'] == dm]
        soilnodes["tag"] = soilnodes.index + 1
        nodetoint = dict(zip(soilnodes['tag'], range(soilnodes.shape[0])))
        cells = soileles[['node1', 'node2', 'node3', 'node4']]
        cells = cells.applymap(nodetoint.get)
        cells = cells.to_numpy(dtype=int)
        soilnodes["Z"] = 0.0
        points = soilnodes[['tag', 'X', 'Y','Z']].to_numpy(dtype=float)
        celltypes = np.ones(cells.shape[0],dtype= int) * pv.CellType.QUAD
        cells = np.insert(cells, 0, 4, axis=1)
        cells[:,1:5] = cells[:,1:5]
        grid = pv.UnstructuredGrid(cells[:,:5], celltypes.tolist(), points[:,1:].tolist())
        pl.set_background('White', top="white")
        # choose matplotlib default colors 
        ss= pl.add_mesh(grid,show_edges=True, cmap="rainbow",style="surface",opacity=1.0)
    pl.show_axes()
    pl.view_xy()
    pl.show()