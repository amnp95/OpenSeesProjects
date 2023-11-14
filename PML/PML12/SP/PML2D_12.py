# %%
import numpy as np 
import pandas as pd 
import Viewmesh2D as vm
import sys
import os
# =============================================================================
# remove all the pervious files
# =============================================================================
os.system('rm boundary*')
os.system('rm *nodes*')
os.system('rm *elements*')
os.system('rm *fixity*')
os.system('rm load*')
os.system('rm *.info')
os.system('rm *.html')

# =======================================================================
# getting odel information from command line
# =======================================================================
Xlen = float(sys.argv[1]);# m
Ylen = float(sys.argv[2]);# m
dx   = float(sys.argv[3]);# m
dy   = float(sys.argv[4]);# m
pmlThickX = float(sys.argv[5]);# m
pmlThickY = float(sys.argv[6]);# m
dxPML = float(sys.argv[7]);# m
dyPML = float(sys.argv[8]);# m
regcores = int(sys.argv[9]);
pmlcore = int(sys.argv[10]);
num_processors = regcores + pmlcore;
eps = 1.e-5;# m

# =======================================================================
# model information
# =======================================================================

# Xlen = 10.;# m
# Ylen = 5. ;# m
# pmlThickX = 3. ;# m
# pmlThickY = 4. ;# m
# dx   = 1. ;# m
# dy   = 1. ;# m
# dxPML = 1. ;# m
# dyPML = .5 ;# m
# eps  = 1.e-5;# m
# regcores = 1;
# pmlcore = 0;
# num_processors = regcores + pmlcore;

# =======================================================================
# creating nodes dataframe
# =======================================================================
# create soil nodes
x = np.arange(-Xlen/2.0,Xlen/2.+eps,dx);
y = np.arange(-Ylen,eps,dy);
Y,X= np.meshgrid(y,x);
nodes = pd.DataFrame({'X':X.flatten(),'Y':Y.flatten()})
nodes["Domain"] = "Soil"


# Adding PML nodes to the nodes dataframe
pmlxlist1 = np.arange(-Xlen/2.0 - pmlThickX, -Xlen/2.0, dxPML)
pmlxlist2 = np.arange(-Xlen/2.0, Xlen/2.0 - eps, dx)
pmlxlist3 = np.arange(Xlen/2.0, Xlen/2.0 + pmlThickX + eps, dxPML)
pmlylist1 = np.arange(-Ylen - pmlThickY, -Ylen, dyPML)
pmlylist2 = np.arange(-Ylen, eps, dy)
pmlxlist = np.concatenate((pmlxlist1, pmlxlist2, pmlxlist3), axis=None)
pmlylist = np.concatenate((pmlylist1, pmlylist2), axis=None)


for x in pmlxlist:
    for y in pmlylist:
        if (x > -Xlen/2.0 +eps) and (x < Xlen/2.0 -eps) and (y > -Ylen +eps) and (y < eps):
            continue
        else:
            nodes_data = {'X':x ,
                          'Y': y,
                        'Domain': 'PML'}
            nodes = pd.concat([nodes, pd.DataFrame([nodes_data])], ignore_index=True)

nodes["Tag"] = nodes.index + 1
# print nodes dataframe to a file
# nodes.to_csv('nodes.csv', index=False)

# =======================================================================
# Creating elements dataframe
# =======================================================================
# create soil elements
elements = pd.DataFrame(columns=[
    'node1', 'node2', 'node3', 'node4'])


# adding soil elements
nx = int(Xlen/dx) 
ny = int(Ylen/dy)
for i in range(nx):
    for j in range(ny):
        node1 = i * (ny + 1) + j + 1
        node2 = node1 + ny + 1
        node3 = node2 + 1;
        node4 = node1 + 1;
        element_data = {'node1': node1, 'node2': node2,
                        'node3': node3, 'node4': node4}
        elements = pd.concat([elements, pd.DataFrame([element_data])], ignore_index=True)

elements["Domain"] = "Soil"

# find maximum nodetag
maxnodetag = nodes['Tag'].max()
# adding PML elements
# create temporary pmlnode datafram from nodes dataframe
pmlnodes = nodes[nodes['Domain'] == 'PML']

for _, row in pmlnodes.iterrows():
    x = row['X']
    y = row['Y']
    dxx = dx
    dyy = dy
    if (x < -Xlen/2.0 -eps) or (x > Xlen/2. -eps):
        dxx = dxPML
    if (y < -Ylen -eps):
        dyy = dyPML
    x1 = x + dxx
    y1 = y + dyy
    n3 = pmlnodes[np.isclose(pmlnodes['X'], x1) & np.isclose(pmlnodes['Y'], y1)].index.values
    if n3.shape[0] == 0:
        # go to the next iteration
        continue

    n4 = pmlnodes[np.isclose(pmlnodes['X'], x ) & np.isclose(pmlnodes['Y'], y1)].index.values
    if n4.shape[0] == 0:
        # go to the next iteration
        continue
    
    
    n2 = pmlnodes[np.isclose(pmlnodes['X'], x1) & np.isclose(pmlnodes['Y'], y)].index.values
    if n2.shape[0] == 0:
        # go to the next iteration
        continue
    
    maxnodetag += 1
    nodes_data = {'X':x + dxx/2.0,
                  'Y':y + dyy/2.0,
                  'Domain': 'PML',
                  'Tag': maxnodetag}
    nodes = pd.concat([nodes, pd.DataFrame([nodes_data])], ignore_index=True)
    
    n1 = row["Tag"]
    n2 = n2[0] + 1
    n3 = n3[0] + 1
    n4 = n4[0] + 1
    n5 = maxnodetag
    element_data = {'node1': n1, 'node2': n2,
                    'node3': n3, 'node4': n4 ,
                    'node5': n5, 'Domain': 'PML'}
    elements = pd.concat([elements, pd.DataFrame([element_data])], ignore_index=True)






# adding tag
elements['tag'] = elements.index + 1
elements['core'] = 0
vm.Domain(nodes, elements)

# =======================================================================
# Creating Soil nodes file
# =======================================================================
# adding status column
nodes['status'] = 0
# create temporary soilnode datafram from nodes dataframe
soilnodes = nodes[nodes['Domain'] == 'Soil']

# iterate over regcores to create a dataframe for each core
for core in range(regcores):
    # filter reg elements
    eles = elements[elements['core'] == core]
    # iterate over elements to node1 to node8 and set their status to 1
    for _, ele in eles.iterrows():
        nodes.loc[ele['node1'] - 1, 'status'] = 1; nodes.loc[ele['node1'] - 1, 'core'] = core
        nodes.loc[ele['node2'] - 1, 'status'] = 1; nodes.loc[ele['node2'] - 1, 'core'] = core
        nodes.loc[ele['node3'] - 1, 'status'] = 1; nodes.loc[ele['node3'] - 1, 'core'] = core
        nodes.loc[ele['node4'] - 1, 'status'] = 1; nodes.loc[ele['node4'] - 1, 'core'] = core

    # filter nodes which are in the core
    soilnodes = nodes[nodes['Domain'] == 'Soil']
    nodes_in_core = soilnodes[soilnodes['status'] == 1]
    
    # write the nodes in a file
    file = open('nodes'+str(core)+'.tcl', 'w')
    for _, node in nodes_in_core.iterrows():
        file.write('node %d %f %f\n' % (node['Tag'], node['X'], node['Y']))
    file.close()

    # reset the status of all nodes to 0
    nodes['status'] = 0


# =======================================================================
# create PML nodes file
# =======================================================================
if num_processors ==1:
    pmlcorelist = [0]
if num_processors>1:
    pmlcorelist = range(regcores, num_processors)
for core in pmlcorelist:
    # filter pml elements
    eles = elements[elements['core'] == core]
    # iterate over elements to node1 to node8 and set their status to 1
    for _, ele in eles.iterrows():
        nodes.loc[ele['node1'] - 1, 'status'] = 1; nodes.loc[ele['node1'] - 1, 'core'] = core
        nodes.loc[ele['node2'] - 1, 'status'] = 1; nodes.loc[ele['node2'] - 1, 'core'] = core
        nodes.loc[ele['node3'] - 1, 'status'] = 1; nodes.loc[ele['node3'] - 1, 'core'] = core
        nodes.loc[ele['node4'] - 1, 'status'] = 1; nodes.loc[ele['node4'] - 1, 'core'] = core
    
    # filter nodes which are in the core
    pmlnodes = nodes[nodes['Domain'] == 'PML']
    nodes_in_core = pmlnodes[pmlnodes['status'] == 1]

    # write the nodes in a file
    file = open('pmlnodes'+str(core)+'.tcl', 'w')
    for _, node in nodes_in_core.iterrows():
        file.write('node %d %f %f\n' % (node['Tag'], node['X'], node['Y']))
    file.close()

    # reset the status of all nodes to 0
    nodes['status'] = 0
# =======================================================================
# create PML center nodes file
# =======================================================================
if num_processors ==1:
    pmlcorelist = [0]
if num_processors>1:
    pmlcorelist = range(regcores, num_processors)

for core in pmlcorelist:
    eles = elements[(elements['core'] == core) & (elements['Domain'] == 'PML')]

    file = open('pmlcenternodes'+str(core)+'.tcl', 'w')
    for _, ele in eles.iterrows():
       file.write('node %d %f %f\n' % (ele['node5'], nodes.loc[ele['node5'] - 1, 'X'], nodes.loc[ele['node5'] - 1, 'Y']))
    file.close()


# =======================================================================
# create Soil elements file
# =======================================================================
# iterate over regcores to create a dataframe for each core
for core in range(regcores):
    eles = elements[(elements['core'] == core) & (elements['Domain'] == 'Soil')]

    file = open('elements'+str(core)+'.tcl', 'w')
    for _, ele in eles.iterrows():
        file.write('eval "element quad %d %d %d %d %d $Soilmaterial"\n' % (ele['tag'], ele['node1'], ele['node2'], ele['node3'], ele['node4']))
    file.close()

# =======================================================================
# create PML elements file
# =======================================================================
for core in pmlcorelist:
    eles = elements[(elements['core'] == core) & (elements['Domain'] == 'PML')]

    file = open('pmlelements'+str(core)+'.tcl', 'w')
    for _, ele in eles.iterrows():
        file.write('eval "element PML12 %d %d %d %d %d %d $PMLmaterial"\n' % (ele['tag'], ele['node1'], ele['node2'], ele['node3'], ele['node4'], ele['node5']))
    file.close()


# =======================================================================
# create boundary file
# =======================================================================
nodes["boundary"] = 0
nodes.loc[(np.isclose(nodes['X'], -Xlen/2., atol=eps) | 
           np.isclose(nodes['X'], Xlen/2., atol=eps)) & 
           (nodes['Y'] >= -Ylen-eps) & (nodes['Y'] <= eps), 'boundary'] = 1

nodes.loc[(np.isclose(nodes['Y'], -Ylen, atol=eps)) &
        (nodes['X'] >= -Xlen/2.-eps) & (nodes['X'] <= +Xlen/2.+eps), 'boundary'] = 1

# seperate the nodes on the with Domain pml and on the bondary
pmlboundarynodes = nodes[(nodes['Domain'] == 'PML') & (nodes['boundary'] == 1)]

# seperate the nodes on the with Domain reg and on the bondary
regboundarynodes = nodes[(nodes['Domain'] == 'Soil') & (nodes['boundary'] == 1)]




# create a map from pmlboundary nodes to regboundarynodes
pmltoreg = dict(zip(pmlboundarynodes.index, regboundarynodes.index))
nodes['status'] = 0
if num_processors ==1:
    pmlcorelist = [0]
if num_processors>1:
    pmlcorelist = range(regcores, num_processors)

for core in pmlcorelist:
    # filter pml elements
    eles = elements[elements['core'] == core]

    # iterate over elements to node1 to node8 and set their status to 1
    for _, ele in eles.iterrows():
        nodes.loc[ele['node1'] - 1, 'status'] = 1
        nodes.loc[ele['node2'] - 1, 'status'] = 1
        nodes.loc[ele['node3'] - 1, 'status'] = 1
        nodes.loc[ele['node4'] - 1, 'status'] = 1

    # extract the list of indexes of nodes which have 1 staus and are on the boundary
    pmlboundarynodes = nodes[(nodes['boundary'] == 1) & (nodes['status'] == 1) &(nodes['Domain'] == 'PML')]
    
    # extract the nodes with the indexes
    # use the pmltoreg map to map the indexes of pmlboundarynodes to regboundarynodes
    regboundarynodes = nodes.loc[pmlboundarynodes.index.map(pmltoreg)]


    # check if the size of the pmlboundarynodes dataframe is equal to the size of the regboundarynodes dataframe and raise an error if not
    if pmlboundarynodes.shape[0] != regboundarynodes.shape[0]:
        raise ValueError('The number of pml boundary nodes is not equal to the number of regular boundary nodes')
    

    # write a file for partioning
    file = open('boundary'+str(core)+'.tcl', 'w')
    for _, (pmlrow, regrow) in enumerate(zip(pmlboundarynodes.iterrows(), regboundarynodes.iterrows())):
    
        # check the coordinates of the nodes are equal or not and raise an error if not
        # using np.isclose because of the floating point precision
        if not (np.isclose(pmlrow[1]['X'], regrow[1]['X'], atol=eps) and np.isclose(pmlrow[1]['Y'], regrow[1]['Y'], atol=eps)):
            raise ValueError('The coordinates of the pml and regular boundary nodes are not equal')
        if pmlcorelist != [0]:
            file.write("node %d %f %f\n" % (regrow[1]['Tag'], regrow[1]['X'], regrow[1]['Y']))
        file.write('equalDOF %d %d 1 2\n' % (regrow[1]['Tag'], pmlrow[1]['Tag']))

    file.close()
    # set the status of the nodes to 0
    nodes['status'] = 0

# =============================================================================
# create loding file
# =============================================================================
# find the nodetags
loadlist = nodes[np.isclose(nodes['X'], 0., atol=eps) & np.isclose(nodes['Y'],  0., atol=eps)]
file = open('load.tcl', 'w')
for _,node in loadlist.iterrows():
    file.write('if {$pid == %d}  {load %d 0 -1}\n' % (0, node['Tag']))


nodetag1 = nodes[np.isclose(nodes['X'], 0.,      atol=eps) &   np.isclose(nodes['Y'],  0.,   atol=eps)]['Tag'].values[0]
nodetag2 = nodes[np.isclose(nodes['X'], Xlen/2., atol=eps) &   np.isclose(nodes['Y'],  0.,   atol=eps)]['Tag'].values[0]
nodetag3 = nodes[np.isclose(nodes['X'], 0,       atol=eps) &   np.isclose(nodes['Y'], -Ylen, atol=eps)]['Tag'].values[0]


file.write ("set recordList {%d %d %d}" % (nodetag1, nodetag2, nodetag3))
file.close()

# %%
