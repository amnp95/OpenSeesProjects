# %%
import numpy as np 
import pandas as pd
import sys

# =============================================================================
# Genreral information about the mesh
# =============================================================================
# getting lx ly lz dx dy dz and pmlthickness from passing arguments to the script
lx = float(sys.argv[1])
ly = float(sys.argv[2])
lz = float(sys.argv[3])
dx = float(sys.argv[4])
dy = float(sys.argv[5])
dz = float(sys.argv[6])
pmlthickness = float(sys.argv[7])

# print the recieved arguments all together
print('lx = %f, ly = %f, lz = %f, dx = %f, dy = %f, dz = %f, pmlthickness = %f' % (lx, ly, lz, dx, dy, dz, pmlthickness))



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


# # print node dataframe to a csv file
# nodes.to_csv('nodes.csv', index=True, header=True)


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
        n1 = pmlnodes[np.isclose(pmlnodes['x'], x) & np.isclose(pmlnodes['y'], y) & np.isclose(pmlnodes['z'], z)].index.values[0]+1
        n2 = pmlnodes[np.isclose(pmlnodes['x'], x1) & np.isclose(pmlnodes['y'], y) & np.isclose(pmlnodes['z'], z)].index.values[0]+1
        n3 = pmlnodes[np.isclose(pmlnodes['x'], x1) & np.isclose(pmlnodes['y'], y1) & np.isclose(pmlnodes['z'], z)].index.values[0]+1
        n4 = pmlnodes[np.isclose(pmlnodes['x'], x) & np.isclose(pmlnodes['y'], y1) & np.isclose(pmlnodes['z'], z)].index.values[0]+1
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


# adding tag
elements['tag'] = elements.index + 1

# print element dataframe to a csv file
# elements.to_csv('elements.csv', index=True, header=True)

# %%
# =============================================================================
# create reg node file 
# =============================================================================
# iterate over the rows of regnodes dataframe to write in file
file = open('nodes.tcl', 'w')
for _, row in nodes[nodes['Domain'] == 'reg'].iterrows():
    file.write('node %d %f %f %f\n' % (row['tag'], row['x'], row['y'], row['z']))
file.close()

# =============================================================================
# create pml node file 
# =============================================================================
# iterate over the rows of pmlnodes dataframe to write in file
file = open('pmlnodes.tcl', 'w')
for _, row in nodes[nodes['Domain'] == 'pml'].iterrows():
    file.write('node %d %f %f %f\n' % (row['tag'], row['x'], row['y'], row['z']))
file.close()


# =============================================================================
# create element file
# =============================================================================
# iterate over the rows of regelements dataframe to write in file
file = open('elements.tcl', 'w')
for _, row in elements[elements['Domain'] == 'reg'].iterrows():
    file.write('element stdBrick %d %d %d %d %d %d %d %d %d $materialTag\n' % (row['tag'], row['node1'], row['node2'], row['node3'], row['node4'], row['node5'], row['node6'], row['node7'], row['node8']))
file.close()


# =============================================================================
# create pml element file
# =============================================================================
# iterate over the rows of pmlelements dataframe to write in file 
file = open('pmlelements.tcl', 'w')
for _, row in elements[elements['Domain'] == 'pml'].iterrows():
    file.write('eval \"element PML %d %d %d %d %d %d %d %d %d $PMLMaterial\"\n' % (row['tag'],row['node1'], row['node2'], row['node3'], row['node4'], row['node5'], row['node6'], row['node7'], row['node8']))
file.close()


# =============================================================================
# create fixity file
# =============================================================================
# iterate over the rows of regnodes dataframe to write in file
file = open('fixity.tcl', 'w')
for _, row in nodes[nodes['Domain'] == 'reg'].iterrows():
    file.write('fix %d 0 1 0\n' % (row['tag']))
file.close()

# iterate over the rows of pmlnodes dataframe to write in file
file = open('pmlfixity.tcl', 'w')
for _, row in nodes[nodes['Domain'] == 'pml'].iterrows():
    file.write('fix %d 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0\n' % (row['tag']))
file.close()
# %%
# =============================================================================
# create boundary file
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

# iterate over pmlboundarynodes and regboundarynodes dataframes and tie the nodes together and write in file
file = open('boundary.tcl', 'w')
for _, (pmlrow, regrow) in enumerate(zip(pmlboundarynodes.iterrows(), regboundarynodes.iterrows())):
    
    # check the coordinates of the nodes are equal or not and raise an error if not
    # using np.isclose because of the floating point precision
    if not (np.isclose(pmlrow[1]['x'], regrow[1]['x'], atol=tolerance) and np.isclose(pmlrow[1]['y'], regrow[1]['y'], atol=tolerance) and np.isclose(pmlrow[1]['z'], regrow[1]['z'], atol=tolerance)):
        raise ValueError('The coordinates of the pml and regular boundary nodes are not equal')
    
    file.write('equalDOF %d %d 1 3\n' % (pmlrow[1]['tag'], regrow[1]['tag']))

file.close()
# %%
# find the node tag of node with corrdinates x = 0. y= 0.5 and z=0. in regnodes dataframe
# using np.isclose because of the floating point precision
nodetag1 = nodes[np.isclose(nodes['x'], 0., atol=tolerance) & np.isclose(nodes['y'],  0.5, atol=tolerance) & np.isclose(nodes['z'], 0, atol=tolerance)]['tag'].values[0]
nodetag2 = nodes[np.isclose(nodes['x'], 0., atol=tolerance) & np.isclose(nodes['y'], -0.5, atol=tolerance) & np.isclose(nodes['z'], 0, atol=tolerance)]['tag'].values[0]
nodetag3 = nodes[np.isclose(nodes['x'], lx/2., atol=tolerance) & np.isclose(nodes['y'], -0.5, atol=tolerance) & np.isclose(nodes['z'], 0, atol=tolerance)]['tag'].values[0]
nodetag4 = nodes[np.isclose(nodes['x'], lx/2., atol=tolerance) & np.isclose(nodes['y'],  0.5, atol=tolerance) & np.isclose(nodes['z'], 0, atol=tolerance)]['tag'].values[0]
nodetag5 = nodes[np.isclose(nodes['x'], 0., atol=tolerance) & np.isclose(nodes['y'],  0.5, atol=tolerance) & np.isclose(nodes['z'], -lz, atol=tolerance)]['tag'].values[0]
nodetag6 = nodes[np.isclose(nodes['x'], 0., atol=tolerance) & np.isclose(nodes['y'],  0.5, atol=tolerance) & np.isclose(nodes['z'], -lz, atol=tolerance)]['tag'].values[0]

# write the node tag in file
file = open('load.tcl', 'w')
file.write('load %d 0 0 -1\n' % (nodetag1))
file.write('load %d 0 0 -1\n' % (nodetag2))
file.write ("set recordList {%d %d %d %d %d %d}" % (nodetag1, nodetag2, nodetag3, nodetag4, nodetag5, nodetag6))
file.close()
