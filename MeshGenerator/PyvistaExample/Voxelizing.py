# %%
from pyvista import examples
import pyvista as pv
import numpy as np
mesh = examples.download_cow()
# %%
mesh.plot(show_edges=True)
# %%
vertices = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0.5, 0.5, -1]])

# mesh faces
faces = np.hstack(
    [
        [4, 0, 1, 2, 3],  # square
        [3, 0, 1, 4],  # triangle
        [3, 1, 2, 4],  # triangle
        [3, 2, 3, 4],  # triangle
        [3, 3, 0, 4],  # triangle
    ]
)

surf = pv.PolyData(vertices, faces)
voxels = pv.voxelize(surf,density=surf.length / 150)
voxels.plot(show_edges=True,scalars=np.arange(voxels.n_points))

# %%
surf.plot(show_edges=True)
# %%
