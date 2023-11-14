# %%
import numpy as np 
import pandas as pd
import viewmesh
import os
import sys
# =============================================================================
# remove all the pervious files
# =============================================================================
os.system('rm boundary*')
os.system('rm *nodes*')
os.system('rm *elements*')
os.system('rm load*')
os.system('rm *.info')
# os.system('rm *.out')
# os.system('rm *.html')
# %%
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
DRM = str(sys.argv[8])



# # # use this for testing
# regcores    = 2
# pmlcores    = 5
# lx          = 10.0
# ly          = 10.0
# lz          = 5.0
# dx          = 1.0
# dy          = 1.0
# dz          = 1.0
# pmlthickness= 2.0


# print the recieved arguments all together
print('lx = %f, ly = %f, lz = %f, dx = %f, dy = %f, dz = %f, pmlthickness = %f' % (lx, ly, lz, dx, dy, dz, pmlthickness, ))



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

pmlxlist = np.arange(xstart - pmlthickness, xend + pmlthickness + eps, dx)
pmlzlist = np.arange(zstart - pmlthickness, zend + eps, dz)
pmlylist = np.arange(ystart - pmlthickness, yend + pmlthickness + eps, dy)


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
            if (x > xstart) and (x < xend) and (y > ystart+eps) and (y < yend-eps) and (z > zstart) and (z < zend + 0.1):
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
# # sort the nodes by x then y then z
# nodes = nodes.sort_values(by=['x', 'y', 'z'])

# # assign indices to the nodes dataframe
# nodes.index = range(0, len(nodes))

# # assing tag to the nodes dataframe
# nodes['tag'] = nodes.index +1

# # # print node dataframe to a csv file
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
# elements.index = range(0, len(elements))

# adding column for element type
elements['Domain'] = 'reg'


# create temporary pmlnode datafram from nodes dataframe
pmlnodes = nodes[nodes['Domain'] == 'pml']
# adding pml elements to the elements dataframe
n1 =1

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
    
    n5 = pmlnodes[np.isclose(pmlnodes['x'], x) & np.isclose(pmlnodes['y'], y) & np.isclose(pmlnodes['z'], z1)].index.values
    if n5.shape[0] == 0:
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

# print element dataframe to a csv file
# elements.to_csv('elements.csv', index=True, header=True)

# %%
# =============================================================================
# create reg node file 
# =============================================================================
# adding status column
nodes['status'] = 0

#iterate over reg nodes to write regnodes file
regnodes = nodes[nodes['Domain'] == 'reg']
file = open('nodes.tcl', 'w')
for _, node in regnodes.iterrows():
    file.write('node %d %f %f %f\n' % (node['tag'], node['x'], node['y'], node['z']))
file.close()
del regnodes
# %%
# =============================================================================
# create pml node file 
# =============================================================================
# iterate over pml nodes to write pmlnodes file
pmlnodes = nodes[nodes['Domain'] == 'pml']
file = open('pmlnodes.tcl', 'w')
for _, node in pmlnodes.iterrows():
    file.write('node %d %f %f %f\n' % (node['tag'], node['x'], node['y'], node['z']))
file.close()
del pmlnodes


# %%
# =============================================================================
# create element file
# =============================================================================
# iterate over reg elements to create element file
regelements = elements[elements['Domain'] == 'reg']
file = open('elements.tcl', 'w')
for _, ele in regelements.iterrows():
    file.write('eval "element stdBrick %d %d %d %d %d %d %d %d %d $materialTag"\n' % (ele['tag'], ele['node1'], ele['node2'], ele['node3'], ele['node4'], ele['node5'], ele['node6'], ele['node7'], ele['node8']))
file.close()
del regelements

# %%
# =============================================================================
# create pml element file
# =============================================================================
# iterate over pml elements to create element file
pmlelements = elements[elements['Domain'] == 'pml']
file = open('pmlelements.tcl', 'w')
for _, ele in pmlelements.iterrows():
    file.write('eval "element PMLVISCOUS %d %d %d %d %d %d %d %d %d $PMLMaterial"\n' % (ele['tag'], ele['node1'], ele['node2'], ele['node3'], ele['node4'], ele['node5'], ele['node6'], ele['node7'], ele['node8']))
file.close()
del pmlelements

# %%
# =============================================================================
# create boundary file
# =============================================================================
# add column says that point is on the boundary of pml or not
# 0 means inside
# 1 means boundary
nodes['boundary'] = 0

tolerance = 1e-6  # Adjust the tolerance based on your needs

# # set the boundary nodes to 1 if x is equal to xstart or xend  z is between zstart and zend
nodes.loc[(np.isclose(nodes['x'], xstart, atol=tolerance) | 
           np.isclose(nodes['x'], xend, atol=tolerance)) & 
           (nodes['y'] >= ystart-eps) & (nodes['y'] <= yend + eps) &
           (nodes['z'] >= zstart-eps) & (nodes['z'] <= zend +eps), 'boundary'] = 1


# set the boundary nodes to 1 if y is equal to ystart or yend  z is between zstart and zend
nodes.loc[(np.isclose(nodes['y'], ystart, atol=tolerance) |
           np.isclose(nodes['y'], yend, atol=tolerance)) &
           (nodes['z'] >= zstart-eps) & (nodes['z'] <= zend + eps) &
           (nodes['x'] >= xstart-eps) & (nodes['x'] <= xend + eps), 'boundary'] = 1

# set the boundary nodes to 1 if z is equal zstart and zend  x is between xstart and xend
nodes.loc[(np.isclose(nodes['z'], zstart, atol=tolerance)) &
        (nodes['x'] >= xstart-eps) & (nodes['x'] <= xend + eps) &
        (nodes['y'] >= ystart-eps) & (nodes['y'] <= yend + eps), 'boundary'] = 1


# seperate the nodes on the with Domain pml and on the bondary
pmlboundarynodes = nodes[(nodes['Domain'] == 'pml') & (nodes['boundary'] == 1)]

# seperate the nodes on the with Domain reg and on the bondary
regboundarynodes = nodes[(nodes['Domain'] == 'reg') & (nodes['boundary'] == 1)]


# checl if the size of the pmlboundarynodes dataframe is equal to the size of the regboundarynodes dataframe and raise an error if not
if pmlboundarynodes.shape[0] != regboundarynodes.shape[0]:
    raise ValueError('The number of pml boundary nodes is not equal to the number of regular boundary nodes')


# create a map from pmlboundary nodes to regboundarynodes
pmltoreg = dict(zip(pmlboundarynodes.index, regboundarynodes.index))



regboundarynodes = nodes.loc[pmlboundarynodes.index.map(pmltoreg)]

# create boundary file
file = open('boundary.tcl', 'w')
for _, (pmlrow, regrow) in enumerate(zip(pmlboundarynodes.iterrows(), regboundarynodes.iterrows())):
        # check the coordinates of the nodes are equal or not and raise an error if not
        # using np.isclose because of the floating point precision
        if not (np.isclose(pmlrow[1]['x'], regrow[1]['x'], atol=tolerance) and np.isclose(pmlrow[1]['y'], regrow[1]['y'], atol=tolerance) and np.isclose(pmlrow[1]['z'], regrow[1]['z'], atol=tolerance)):
            raise ValueError('The coordinates of the pml and regular boundary nodes are not equal')
        # file.write('equalDOF %d %d 1 2 3\n' % (regrow[1]['tag'], pmlrow[1]['tag']))
        file.write('equalDOF %d %d 1 2 3\n' % (pmlrow[1]['tag'], regrow[1]['tag']))
        # file.write("rigidLink bar %d %d\n" % (regrow[1]['tag'], pmlrow[1]['tag']))
        # file.write("rigidLink bar %d %d\n" % (pmlrow[1]['tag'], regrow[1]['tag']))


 
# %%
# =============================================================================
# create loding file
# =============================================================================
# find the nodetags
# loadlist = nodes[(np.abs(nodes['x']) < 1+eps)  & (np.abs(nodes['y']) < 1+eps) & np.isclose(nodes['z'], 0, atol=tolerance)]
loadlist = nodes[np.isclose(nodes['x'], 0., atol=tolerance) & np.isclose(nodes['y'],  0., atol=tolerance) & np.isclose(nodes['z'], 0, atol=tolerance)]
nodetag1 = nodes[np.isclose(nodes['x'], 0., atol=tolerance) & np.isclose(nodes['y'],  0., atol=tolerance) & np.isclose(nodes['z'], 0, atol=tolerance)]['tag'].values[0]
# nodetag2 = nodes[np.isclose(nodes['x'], lx/2., atol=tolerance) & np.isclose(nodes['y'],  0., atol=tolerance) & np.isclose(nodes['z'], 0, atol=tolerance)]['tag'].values[0]
# nodetag3 = nodes[np.isclose(nodes['x'], 0, atol=tolerance) & np.isclose(nodes['y'],  ly/2., atol=tolerance) & np.isclose(nodes['z'], 0, atol=tolerance)]['tag'].values[0]
# nodetag4 = nodes[np.isclose(nodes['x'], lx/2., atol=tolerance) & np.isclose(nodes['y'],  ly/2., atol=tolerance) & np.isclose(nodes['z'], 0, atol=tolerance)]['tag'].values[0]
# nodetag5 = nodes[np.isclose(nodes['x'], 0., atol=tolerance) & np.isclose(nodes['y'],  0., atol=tolerance) & np.isclose(nodes['z'], -lz+1, atol=tolerance)]['tag'].values[0]
# nodetag6 = nodes[np.isclose(nodes['x'], lx/2., atol=tolerance) & np.isclose(nodes['y'],  0., atol=tolerance) & np.isclose(nodes['z'], -lz+1, atol=tolerance)]['tag'].values[0]
# nodetag7 = nodes[np.isclose(nodes['x'], 0, atol=tolerance) & np.isclose(nodes['y'],  ly/2., atol=tolerance) & np.isclose(nodes['z'], -lz+1, atol=tolerance)]['tag'].values[0]
# nodetag8 = nodes[np.isclose(nodes['x'], lx/2., atol=tolerance) & np.isclose(nodes['y'],  ly/2., atol=tolerance) & np.isclose(nodes['z'], -lz+1, atol=tolerance)]['tag'].values[0]


# write the node tag in file
file = open('load.tcl', 'w')
for _,node in loadlist.iterrows():
    file.write('load %d 0 0 -1\n' % (node['tag']))
file.close()

file = open("record.tcl", "w")
file.write ("set recordList {%d}" % (nodetag1))
file.close()



# %%
# =============================================================================
# Create DRM node tag file
# =============================================================================
if DRM == "YES":
    # read the coords
    coords = np.loadtxt("/home/amnp95/Projects/Github/OpenSeesProjects/Shakermaker/SP/DRMLOAD/coords.csv", delimiter=",")
    file = open("DRMLOAD/nodeTags.csv", "w")
    file2 = open("load.tcl","w")
    file2.write("set loadList {")
    print(coords.shape)
    for i in range(coords.shape[0]):
        tags = nodes.loc[(np.isclose(nodes['x'], coords[i,0], atol=tolerance) & np.isclose(nodes['y'],  coords[i,1], atol=tolerance) & np.isclose(nodes['z'], coords[i,2], atol=tolerance))]["tag"].values[0]
        file.write("%d\n" % (tags))
        file2.write("%d " % (tags))
    file2.write("}")
    file2.close()
    file.close()




# %%
# =============================================================================
# plot the mesh with deifferent cores
# =============================================================================
# read the tags from the tags.txt file
tags = np.loadtxt('tags.txt', dtype=int)
# set the core of the elements with tag in tags list to 1
viewmesh.view(nodes.copy(), elements.copy(), tags)
# viewmesh.cores(nodes.copy(), elements.copy(),  view="pml")
# viewmesh.cores(nodes.copy(), elements.copy(),  view="reg")
viewmesh.cores(nodes.copy(), elements.copy(),  view="pml")
# viewmesh.cores(nodes.copy(), elements.copy(),  view="all")


# %%
