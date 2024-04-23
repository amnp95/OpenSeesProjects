# # =======================================================================================================
# # add building 
# # =======================================================================================================
# # %%
# pl = pv.Plotter(off_screen=False)
# pl.set_background('white')
# xfloor  = np.arange(-8, 8, 2.5, dtype=np.float32)
# yfloor  = np.arange(-8, 8, 2.5, dtype=np.float32)
# zfloor  = np.arange(0, 1.1, 1., dtype=np.float32)
# x,y,z   = np.meshgrid(xfloor,yfloor,zfloor,indexing='ij')
# floor   = pv.StructuredGrid(x,y,z)
# floors  = floor.copy()
# floors.points[:,2] += 5

# for i in range(1,8):
#     f = floor.copy()
#     f.points[:,2] += 5  + 5 * i
    
#     floors = floors.merge(f)

# # adding columns
# X = np.arange(-7, 7, 4, dtype=np.float32)
# Y = np.arange(-7, 7, 4, dtype=np.float32)
# height = 40
# count = 0
# for x in X:
#     for y in Y:
#         clindeer = pv.Cylinder(center=(x,y,height/2.-2.5), direction=(0,0,1), radius=0.25, height=height+5)
#         if count == 0:
#             columns = clindeer
#         else:
#             columns = columns.merge(clindeer)
#         count += 1



# for detail in [columns, floors]:
#     detail.points[:,0] += 30
#     detail.points[:,1] += 50
#     detail.points[:,2] += 20

# # pl.add_mesh(floors, show_edges=False, style="surface", opacity=1.0)
# # pl.add_mesh(columns, show_edges=False, style="surface", opacity=1.0, color="purple")
# mesh2 = mesh.copy()
# mesh2.clear_cell_data()
# mesh2.clear_point_data()
# matplotlib_defaultcolors.append("navy")
# matplotlib_defaultcolors.append("khaki")
# matplotlib_defaultcolors.append("olive")
# # pl.add_mesh(regular,scalars="partitioned", show_edges=True,opacity=1.0,cmap=matplotlib_defaultcolors[:reg_num_cores])
# pl.add_mesh(PML,scalars="partitioned",show_edges=True,opacity=1.0,color=matplotlib_defaultcolors[1],cmap=matplotlib_defaultcolors[reg_num_cores+DRM_num_cores:reg_num_cores+DRM_num_cores+PML_num_cores])
# # pl.add_mesh(DRM,scalars="partitioned", show_edges=True,opacity=1.0,color=matplotlib_defaultcolors[2],cmap=matplotlib_defaultcolors[reg_num_cores:reg_num_cores+DRM_num_cores])
# pl.remove_scalar_bar()
# pl.enable_anti_aliasing('ssaa')

# # pl.add_mesh(columns, show_edges=False, style="surface", opacity=1.0, color="red")
# pl.show()