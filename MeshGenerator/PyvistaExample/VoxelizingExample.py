import numpy as np

import pyvista as pv

# sphinx_gallery_thumbnail_number = 2
from pyvista import examples

# Load a surface to voxelize
surface = examples.download_foot_bones()
surface.plot()
# surface.plot()
voxels = pv.voxelize(surface, density=surface.length / 200)

p = pv.Plotter()
p.add_mesh(voxels, color=True, show_edges=True, opacity=1.0)
# p.add_mesh(surface, color="lightblue", opacity=0.5)
p.show()
###############################################################################
# We could even add a scalar field to that new voxel model in case we
# wanted to create grids for modelling. In this case, let's add a scalar field
# for bone density noting:
voxels["density"] = np.full(voxels.n_cells, 3.65)  # g/cc
voxels

###############################################################################
voxels.plot(scalars="density")


###############################################################################
# A constant scalar field is kind of boring, so let's get a little fancier by
# added a scalar field that varies by the distance from the bounding surface.
voxels.compute_implicit_distance(surface, inplace=True)
voxels

###############################################################################
contours = voxels.contour(6, scalars="implicit_distance")

p = pv.Plotter()
# p.add_mesh(voxels, opacity=0.25, scalars="implicit_distance")
p.add_mesh(contours, opacity=0.5, scalars="implicit_distance")
p.show()