import pyvista as pv
import numpy as np

# create a function for view mesh in pyvista

def viewcores(nodes, elements):
    # =============================================================================
    # plot the mesh with deifferent cores
    # =============================================================================
    # set nodes status to 0
    nodes['status'] = 0


    # initialize plotter and its color cycler
    pl = pv.Plotter()
    pl.set_color_cycler("default")

    # find number of cores form the core column in elements
    cores = elements['core'].max() + 1

    for core in range(cores):
        
        # filter elements
        eles = elements[elements['core'] == core]

        # create cell array and point array
        # cell should be integer and point should be float


        # iterate over elements to node1 to node8 and set their status to 1
        for _, ele in eles.iterrows():
            nodes.loc[ele['node1'] - 1, 'status'] = 1
            nodes.loc[ele['node2'] - 1, 'status'] = 1
            nodes.loc[ele['node3'] - 1, 'status'] = 1
            nodes.loc[ele['node4'] - 1, 'status'] = 1
            nodes.loc[ele['node5'] - 1, 'status'] = 1
            nodes.loc[ele['node6'] - 1, 'status'] = 1
            nodes.loc[ele['node7'] - 1, 'status'] = 1
            nodes.loc[ele['node8'] - 1, 'status'] = 1
        
        # filter nodes
        nodes_in_core = nodes[nodes['status'] == 1]

        # create map between the node tag and the integer between 0 and the number of nodes in the core
        nodetoint = dict(zip(nodes_in_core['tag'], range(nodes_in_core.shape[0])))

        # map the nodes in the cells to the integer
        cells = eles[['node1', 'node2', 'node3', 'node4', 'node5', 'node6', 'node7', 'node8']]

        # map the nodes to the integer
        cells = cells.applymap(nodetoint.get)
    
        # convert the cells to numpy array
        cells = cells.to_numpy(dtype=int)

        # create point array
        points = nodes_in_core[['tag', 'x', 'y', 'z']].to_numpy(dtype=float)
        celltypes = np.ones(cells.shape[0],dtype= int) * pv.CellType.HEXAHEDRON
        # add cloumn of eight to cells at the begining
        cells = np.insert(cells, 0, 8, axis=1)

        cells[:,1:9] = cells[:,1:9]
        grid = pv.UnstructuredGrid(cells[:,:9], celltypes.tolist(), points[:,1:].tolist())
        pl.set_background('White', top="white")
        # choose matplotlib default colors 
        
        ss= pl.add_mesh(grid,show_edges=True, cmap="rainbow",style="surface",opacity=1.0)
        # save the screenshot
    pl.show_axes()
    pl.show()