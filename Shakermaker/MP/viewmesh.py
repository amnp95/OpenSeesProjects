import pyvista as pv
import numpy as np


def cores (nodes, elements, regcores, pmlcores, view="all") :
    '''
    View sparsed the mesh of the 3D box model.

    Parameters
    ----------
    nodes : pandas.DataFrame 
        The nodes of the mesh.

    elements : pandas.DataFrame
        The elements of the mesh.

    regcores : int
        The number of regular cores.

    pmlcores : int
        The number of PML cores.

    view : str, optional
        The view of the mesh. The default is "all".
        can be regular or pml or an integer between 0 and the number of cores.
    '''
    # initialize plotter and its color cycler
    pl = pv.Plotter()
    pl.set_color_cycler("default")

    # check if view is integer
    if isinstance(view, int):
        start = view
        end = view + 1
    elif view == "regular":
        start = 0
        end = regcores
    elif view == "pml":
        start = regcores
        end = regcores + pmlcores
    elif view == "all":
        start = 0
        end = regcores + pmlcores


    for core in range(start, end):
    
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

        # print(nodes_in_core)
        # print(cells)


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
        # pl.export_html(
    #     'cores.html', backend='panel'
    # )  
    pl.show_axes()
    pl.show()
 