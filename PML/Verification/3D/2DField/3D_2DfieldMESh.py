# %%
import numpy as np 
import pandas as pd


# =============================================================================
# Genreral information about the mesh
# =============================================================================
xstart = -3.
ystart = -0.5
zstart = -3.

xend   = 3.
yend   = 0.5
zend   = 0.

dx     = 1.0
dy     = 1.0
dz     = 1.0
eps   = 1e-6

x = np.arange(xstart,xend + eps, dx)
y = np.arange(ystart,yend + eps, dy)
z = np.arange(zstart,zend + eps, dz)

nx = len(x)-1
ny = len(y)-1
nz = len(z)-1

pmlthickness = 2.0

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


# print node dataframe to a csv file
nodes.to_csv('nodes.csv', index=True, header=True)


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

# adding column for element type
elements['Domain'] = 'reg'



# create temporary pmlnode datafram from nodes dataframe
pmlnodes = nodes[nodes['Domain'] == 'pml']
# adding pml elements to the elements dataframe



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
    

    if pmlnodes[np.isclose(pmlnodes['x'], x1) & np.isclose(pmlnodes['y'], y1) & np.isclose(pmlnodes['z'], z1)].index.values.shape[0] == 1:
        
        # searching for 8 nodes of the element with 
        n1 = pmlnodes[np.isclose(pmlnodes['x'], x) & np.isclose(pmlnodes['y'], y) & np.isclose(pmlnodes['z'], z)].index.values[0]
        n2 = pmlnodes[np.isclose(pmlnodes['x'], x1) & np.isclose(pmlnodes['y'], y) & np.isclose(pmlnodes['z'], z)].index.values[0]
        n3 = pmlnodes[np.isclose(pmlnodes['x'], x1) & np.isclose(pmlnodes['y'], y1) & np.isclose(pmlnodes['z'], z)].index.values[0]
        n4 = pmlnodes[np.isclose(pmlnodes['x'], x) & np.isclose(pmlnodes['y'], y1) & np.isclose(pmlnodes['z'], z)].index.values[0]
        n5 = n1 + 1
        n6 = n2 + 1
        n7 = n3 + 1
        n8 = n4 + 1


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




# print element dataframe to a csv file
# elements.to_csv('elements.csv', index=True, header=True)

# %%
# =============================================================================
# create reg node file 
# =============================================================================
regnodes = nodes[nodes['Domain'] == 'reg']

# change the reg to "node"
regnodes['Domain'] = 'node'

# relocate node column to the first column
cols = regnodes.columns.tolist()
cols = cols[-2:] + cols[:-2]
regnodes = regnodes[cols]


# print regnodes dataframe in a "tag xcoord ycoord zcoord" format
regnodes.to_csv('nodes.tcl', sep=' ',index=False, header=False)


# =============================================================================
# create pml node file 
# =============================================================================
regnodes = nodes[nodes['Domain'] == 'pml']

# change the reg to "node"
regnodes['Domain'] = 'node'

# relocate node column to the first column
cols = regnodes.columns.tolist()
cols = cols[-2:] + cols[:-2]
regnodes = regnodes[cols]


# print regnodes dataframe in a "tag xcoord ycoord zcoord" format
regnodes.to_csv('pmlnodes.tcl', sep=' ',index=False, header=False)

# unaasigned regnodes to free up memory
del regnodes

# =============================================================================
# create element file
# =============================================================================
regelements = elements[elements['Domain'] == 'reg']


# change the reg to "element "
regelements['Domain'] = 'element'

# adding type column
regelements['type'] = 'stdBrick'

# adding tag
elements['tag'] = elements.index + 1

#relocate node column to the first column
cols = regelements.columns.tolist()
cols = cols[-2:] + cols[:-2]
regelements = regelements[cols]

# print regelements dataframe in a "tag xcoord ycoord zcoord" format
regelements.to_csv('elements.tcl', sep=' ',index=False, header=False)


# unaasigned regelements to free up memory
del regelements

# =============================================================================
# create pml element file
# =============================================================================
pmlelements = elements[elements['Domain'] == 'pml']

# change the reg to "element "
pmlelements['Domain'] = 'element PML'

# relocate node column to the first column
cols = pmlelements.columns.tolist()
cols = cols[-2:] + cols[:-2]
pmlelements = pmlelements[cols]


# print pmlelements dataframe in a "tag xcoord ycoord zcoord" format
pmlelements.to_csv('pmlelements.tcl', sep=' ',index=False, header=False)


# unaasigned pmlelements to free up memory
del pmlelements

# %%

# add column says that point is on the boundary or not
# 0 means inside
# 1 means boundary


nodes['boundary'] = 0


# tolerance = 1e-6  # Adjust the tolerance based on your needs
# nodes.loc[
#     (np.isclose(nodes['x'], xstart, atol=tolerance) |
#      np.isclose(nodes['x'], xend, atol=tolerance) |
#      np.isclose(nodes['y'], ystart, atol=tolerance) |
#      np.isclose(nodes['y'], yend, atol=tolerance) |
#      np.isclose(nodes['z'], zstart, atol=tolerance) |
#      np.isclose(nodes['z'], zend, atol=tolerance)),
#     'boundary'
# ] = 1
# %%
