# %%
import numpy as np 
import pandas as pd
import pyvista as pv
import viewmesh as vm
import os
import sys
# =============================================================================
# remove all the pervious files
# =============================================================================
# %%
os.system('rm boundary*')
os.system('rm *nodes*')
os.system('rm *elements*')
os.system('rm *fixity*')
os.system('rm load*')
os.system('rm *.info')
os.system('rm *.html')
# %%
# =============================================================================
# Genreral information about the mesh
# =============================================================================
# getting lx ly lz dx dy dz and pmlthickness from passing arguments to the script
# cores        = int(sys.argv[1])
# lx           = float(sys.argv[2])
# ly           = float(sys.argv[3])
# lz           = float(sys.argv[4])
# dx           = float(sys.argv[5])
# dy           = float(sys.argv[6])
# dz           = float(sys.argv[7])
# pmlthickness = float(sys.argv[8])

# use this for testing
lx          = 20.0
ly          = 1.0
lz          = 10.0
dx          = 1.0
dy          = 1.0
dz          = 1.0
pmlthickness= 2.0
cores       = 2

# print the recieved arguments all together
print('cores = %d, lx = %f, ly = %f, lz = %f, dx = %f, dy = %f, dz = %f, pmlthickness = %f' % (cores, lx, ly, lz, dx, dy, dz, pmlthickness))

xstart = -lx/2
ystart = -ly/2
zstart = -lz

xend   = lx/2
yend   = ly/2
zend   = 0.0

eps   = 1e-6

x = np.arange(xstart,xend + eps, dx)
y = np.arange(ystart,yend + eps, dy)
z = np.arange(zstart,zend + eps, dz)

nx = len(x)-1
ny = len(y)-1
nz = len(z)-1

pmlxlist = np.arange(xstart - pmlthickness, xend +pmlthickness + eps, dx)
pmlzlist = np.arange(zstart - pmlthickness, zend + eps, dz)
pmlylist = np.arange(ystart , yend + eps, dy)

# %%
# =============================================================================
# create nodes dataframe
# =============================================================================
# create a meshgrid

X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

# create a dataframe
nodes = pd.DataFrame({'x':X.flatten(),'y':Y.flatten(),'z':Z.flatten()})
X = Y = Z = None


# adding column for the domain
nodes['Domain'] = 'reg'



# Adding PML nodes to the nodes dataframe
for x in pmlxlist:
    for y in pmlylist:
        for z in pmlzlist:
            if (x > xstart) and (x < xend) and (y > ystart-eps) and (y < yend+eps) and (z > zstart) and (z < zend + 0.1):
                continue
            else:
                nodes_data = {'x':x ,
                              'y': y,
                              'z': z,
                              'Domain': 'pml'}
                
                nodes = pd.concat([nodes, pd.DataFrame([nodes_data])], ignore_index=True)


# adding tag 
nodes['tag'] = nodes.index +1
nodes ['core'] = 0


# %%
# =============================================================================
# Create a dataframe for elements
# =============================================================================
# Create an empty dataframe
# 8 nodes for each element

elements = pd.DataFrame(columns=[
    'node1', 'node2', 'node3', 'node4', 'node5', 'node6', 'node7', 'node8'])


# Add nodes to each element
for x in range(nx):
    for y in range(ny):
        for z in range(nz):
            node1 = int(x * (ny + 1) * (nz + 1) + y * (nz + 1) + z + 1)
            node2 = int((x + 1) * (ny + 1) * (nz + 1) + y * (nz + 1) + z + 1)
            node3 = int((x + 1) * (ny + 1) * (nz + 1) + (y + 1) * (nz + 1) + z + 1)
            node4 = int(x * (ny + 1) * (nz + 1) + (y + 1) * (nz + 1) + z + 1)
            node5 = node1 + 1
            node6 = node2 + 1
            node7 = node3 + 1
            node8 = node4 + 1
            

            element_data = {
                'node1': node1,
                'node2': node2,
                'node3': node3,
                'node4': node4,
                'node5': node5,
                'node6': node6,
                'node7': node7,
                'node8': node8,
            }

            elements = pd.concat([elements, pd.DataFrame([element_data])], ignore_index=True)


# Assign indices to the elements dataframe
elements.index = range(1, len(elements) + 1)
# elements.index = range(0, len(elements))

# adding column for element type
elements['Domain'] = 'reg'


# create temporary pmlnode datafram from nodes dataframe
pmlnodes = nodes[nodes['Domain'] == 'pml']
# adding pml elements to the elements dataframe
n1 = 1

# Iterate over the rows of pmlnodes dataframe
for _, row in pmlnodes.iterrows():
    x = row['x']
    y = row['y']
    z = row['z']
    x1 = x + dx
    y1 = y + dy
    z1 = z + dz

    # check that if a node exist in the pmlnodes dataframe with x1 y1 z1 coordinates
    # if pmlnodes[np.isclose(pmlnodes['x'], x1) & np.isclose(pmlnodes['y'], y1) & np.isclose(pmlnodes['z'], z1)] is not empty
    
    # using np.isclose because of the floating point precision
 
    n7 = pmlnodes[np.isclose(pmlnodes['x'], x1) & np.isclose(pmlnodes['y'], y1) & np.isclose(pmlnodes['z'], z1)].index.values

    if n7.shape[0] == 0:
        # go to the next iteration
        continue
    
    n6 = pmlnodes[np.isclose(pmlnodes['x'], x1) & np.isclose(pmlnodes['y'], y) & np.isclose(pmlnodes['z'], z1)].index.values
    if n6.shape[0] == 0:
        # go to the next iteration
        continue

    n8 = pmlnodes[np.isclose(pmlnodes['x'], x) & np.isclose(pmlnodes['y'],  y1) & np.isclose(pmlnodes['z'], z1)].index.values
    if n8.shape[0] == 0:
        # go to the next iteration
        continue
    
    n1 = row['tag']
    n8 = n8[0] + 1
    n6 = n6[0] + 1
    n7 = n7[0] + 1
    n5 = n1 + 1
    n2 = n6 - 1
    n3 = n7 - 1
    n4 = n8 - 1


    element_data = {
        'node1': n1,
        'node2': n2,
        'node3': n3,
        'node4': n4,
        'node5': n5,
        'node6': n6,
        'node7': n7,
        'node8': n8,
        'Domain': 'pml'
    }
    elements = pd.concat([elements, pd.DataFrame([element_data])], ignore_index=True)


# adding tag
elements['tag'] = elements.index + 1

# adding core column
elements['core'] = 0

# %%
# =============================================================================
# sparse elements to different cores based on geometry
# =============================================================================

for _, row in elements.iterrows():
    node1 = row['node1']
    if nodes.loc[node1 -1, 'x'] < 0:
        elements.loc[_, 'core'] = 1

 
vm.viewcores(nodes.copy(), elements.copy())
 

# %%
# =============================================================================
# create node file 
# =============================================================================
# adding status column
nodes['status'] = 0

# iterate over regcores to create a dataframe for each core
for core in range(cores):
    # filter reg elements
    eles = elements[elements['core'] == core]
    # iterate over elements to node1 to node8 and set their status to 1
    for _, ele in eles.iterrows():
        nodes.loc[ele['node1'] - 1, 'status'] = 1; nodes.loc[ele['node1'] - 1, 'core'] = core
        nodes.loc[ele['node2'] - 1, 'status'] = 1; nodes.loc[ele['node2'] - 1, 'core'] = core
        nodes.loc[ele['node3'] - 1, 'status'] = 1; nodes.loc[ele['node3'] - 1, 'core'] = core
        nodes.loc[ele['node4'] - 1, 'status'] = 1; nodes.loc[ele['node4'] - 1, 'core'] = core
        nodes.loc[ele['node5'] - 1, 'status'] = 1; nodes.loc[ele['node5'] - 1, 'core'] = core
        nodes.loc[ele['node6'] - 1, 'status'] = 1; nodes.loc[ele['node6'] - 1, 'core'] = core
        nodes.loc[ele['node7'] - 1, 'status'] = 1; nodes.loc[ele['node7'] - 1, 'core'] = core
        nodes.loc[ele['node8'] - 1, 'status'] = 1; nodes.loc[ele['node8'] - 1, 'core'] = core


    # filter nodes which are in the core
    nodes_in_core = nodes[nodes['status'] == 1]

    file = open('nodes'+str(core)+'.tcl', 'w')
    
    # write the reg nodes in a file
    reg_nodes = nodes_in_core[nodes_in_core['Domain'] == 'reg']
    file.write("model BasicBuilder -ndm 3 -ndf 3\n")
    for _, node in reg_nodes.iterrows():
        file.write('node %d %f %f %f\n' % (node['tag'], node['x'], node['y'], node['z']))
        file.write('fix %d 0 1 0\n' % (node['tag']))

    pml_nodes = nodes_in_core[nodes_in_core['Domain'] == 'pml']
    # write the nodes in a file

    # reset the status of all nodes to 0
    file.write("model BasicBuilder -ndm 3 -ndf 18\n")
    for _, node in pml_nodes.iterrows():
        file.write('node %d %f %f %f\n' % (node['tag'], node['x'], node['y'], node['z']))
        file.write('fix %d 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0\n' % (node['tag']))
    file.close()

    nodes['status'] = 0

# delte the status column
del nodes['status']

# delete eles and nodes_in_core
# check if the eles dataframe exists
if 'eles' in locals():
    del eles, nodes_in_core


# %%
# =============================================================================
# create element file
# =============================================================================
# iterate over each core to create element file
for core in range(cores):
    # filter reg elements
    eles = elements[elements['core'] == core]

    regeles = eles[eles['Domain'] == 'reg']
    # write a file for partioning
    file = open('elements'+str(core)+'.tcl', 'w')
    
    # iterate over the rows of reg elements dataframe to write in file
    for _, row in regeles.iterrows():
        file.write('element stdBrick %d %d %d %d %d %d %d %d %d $materialTag\n' % (row['tag'], row['node1'], row['node2'], row['node3'], row['node4'], row['node5'], row['node6'], row['node7'], row['node8']))
    
    # filter pml elements
    pml_eles = eles[eles['Domain'] == 'pml']
    
    # iterate over the rows of pml elements dataframe to write in file
    for _, row in pml_eles.iterrows():
        file.write('eval "element PML %d %d %d %d %d %d %d %d %d $PMLMaterial"\n' % (row['tag'],row['node1'], row['node2'], row['node3'], row['node4'], row['node5'], row['node6'], row['node7'], row['node8']))
    file.close()



# %%
# =============================================================================
# create mp constraint file
# =============================================================================
# add column says that point is on the boundary of pml or not
# 0 means inside
# 1 means boundary
nodes['boundary'] = 0

tolerance = 1e-6  # Adjust the tolerance based on your needs

# set the boundary nodes to 1 if x is equal to xstart or xend  z is between zstart and zend
nodes.loc[(np.isclose(nodes['x'], xstart, atol=tolerance) | np.isclose(nodes['x'], xend, atol=tolerance)) & (nodes['z'] >= zstart-eps) & (nodes['z'] <= zend +eps), 'boundary'] = 1

# set the boundary nodes to 1 if z is equal zstart and zend  x is between xstart and xend
nodes.loc[(np.isclose(nodes['z'], zstart, atol=tolerance)) & (nodes['x'] >= xstart-eps) & (nodes['x'] <= xend +eps), 'boundary'] = 1



# seperate the nodes on the with Domain pml and on the bondary
pmlboundarynodes = nodes[(nodes['Domain'] == 'pml') & (nodes['boundary'] == 1)]

# seperate the nodes on the with Domain reg and on the bondary
regboundarynodes = nodes[(nodes['Domain'] == 'reg') & (nodes['boundary'] == 1)]

# checl if the size of the pmlboundarynodes dataframe is equal to the size of the regboundarynodes dataframe and raise an error if not
if pmlboundarynodes.shape[0] != regboundarynodes.shape[0]:
    raise ValueError('The number of pml boundary nodes is not equal to the number of regular boundary nodes')


# create a map from pmlboundary nodes to regboundarynodes
pmltoreg = dict(zip(pmlboundarynodes.index, regboundarynodes.index))


# add status column to nodes dataframe
nodes['status'] = 0

# iterate over each core to create boundary file
for core in range(cores):
    # filter pml elements
    eles = elements[elements['core'] == core]

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

    # extract the list of indexes of nodes which have 1 status and are on the boundary
    pmlboundarynodes = nodes[(nodes['boundary'] == 1) & (nodes['status'] == 1) & (nodes['Domain'] == 'pml') ]
    # extract the nodes with the indexes
    # use the pmltoreg map to map the indexes of pmlboundarynodes to regboundarynodes
    regboundarynodes = nodes.loc[pmlboundarynodes.index.map(pmltoreg)]


    # check if the size of the pmlboundarynodes dataframe is equal to the size of the regboundarynodes dataframe and raise an error if not
    if pmlboundarynodes.shape[0] != regboundarynodes.shape[0]:
        raise ValueError('The number of pml boundary nodes is not equal to the number of regular boundary nodes')
    

    # write a file for partioning
    file = open('mpconstraint'+str(core)+'.tcl', 'w')
    for _, (pmlrow, regrow) in enumerate(zip(pmlboundarynodes.iterrows(), regboundarynodes.iterrows())):
    
        # check the coordinates of the nodes are equal or not and raise an error if not
        # using np.isclose because of the floating point precision
        if not (np.isclose(pmlrow[1]['x'], regrow[1]['x'], atol=tolerance) and np.isclose(pmlrow[1]['y'], regrow[1]['y'], atol=tolerance) and np.isclose(pmlrow[1]['z'], regrow[1]['z'], atol=tolerance)):
            raise ValueError('The coordinates of the pml and regular boundary nodes are not equal')
        file.write('equalDOF %d %d 1 3\n' % (pmlrow[1]['tag'], regrow[1]['tag']))

    file.close()
    # set the status of the nodes to 0
    nodes['status'] = 0
# =============================================================================
# create loding file
# =============================================================================
# find the nodetags
nodetag1 = nodes[np.isclose(nodes['x'], 0., atol=tolerance) & np.isclose(nodes['y'],  0.5, atol=tolerance) & np.isclose(nodes['z'], 0, atol=tolerance)]
nodetag2 = nodes[np.isclose(nodes['x'], 0., atol=tolerance) & np.isclose(nodes['y'], -0.5, atol=tolerance) & np.isclose(nodes['z'], 0, atol=tolerance)]
nodetag3 = nodes[np.isclose(nodes['x'], lx/2., atol=tolerance) & np.isclose(nodes['y'], -0.5, atol=tolerance) & np.isclose(nodes['z'], 0, atol=tolerance)]['tag'].values[0]
nodetag4 = nodes[np.isclose(nodes['x'], lx/2., atol=tolerance) & np.isclose(nodes['y'],  0.5, atol=tolerance) & np.isclose(nodes['z'], 0, atol=tolerance)]['tag'].values[0]
nodetag5 = nodes[np.isclose(nodes['x'], 0., atol=tolerance) & np.isclose(nodes['y'],  0.5, atol=tolerance) & np.isclose(nodes['z'], -lz, atol=tolerance)]['tag'].values[0]
nodetag6 = nodes[np.isclose(nodes['x'], 0., atol=tolerance) & np.isclose(nodes['y'],  0.5, atol=tolerance) & np.isclose(nodes['z'], -lz, atol=tolerance)]['tag'].values[0]


# write the node tag in file
file = open('load.tcl', 'w')
file.write('if {$pid == %d}  {load %d 0 0 -1}\n' % (nodetag1['core'].values[0], nodetag1['tag'].values[0]))
file.write('if {$pid == %d}  {load %d 0 0 -1}\n' % (nodetag2['core'].values[0], nodetag2['tag'].values[0]))
file.write ("set recordList {%d %d %d %d %d %d}" % (nodetag1['tag'].values[0], nodetag2['tag'].values[0], nodetag3, nodetag4, nodetag5, nodetag6))
file.close()










# # %%
# # view the mesh with pyvista
# # create cell array and point array 
# # cell should be integer and point should be float
# cells = elements[['tag', 'node1', 'node2', 'node3', 'node4', 'node5', 'node6', 'node7', 'node8', 'core']].to_numpy(dtype=int)
# points = nodes[['tag', 'x', 'y', 'z']].to_numpy(dtype=float)


# celltypes = np.ones(cells.shape[0],dtype= int) * pv.CellType.HEXAHEDRON
# cells[:,0] = 8
# cells[:,1:9] = cells[:,1:9] -1
# grid = pv.UnstructuredGrid(cells[:,:9], celltypes.tolist(), points[:,1:].tolist())
# grid.cell_data["processors"] = cells [:,-1]
# pl = pv.Plotter()
# pl.set_background('White', top="white")
# ss= pl.add_mesh(grid,show_edges=True, cmap="rainbow",style="surface",opacity=1.0)
# pl.show()

