import numpy as np

import pyvista as pv
from pyvista import CellType
lin_pts = np.array(
    [
        [-1, -1, -1],  # point 0
        [1, -1, -1],  # point 1
        [1, 1, -1],  # point 2
        [-1, 1, -1],  # point 3
        [-1, -1, 1],  # point 4
        [1, -1, 1],  # point 5
        [1, 1, 1],  # point 6
        [-1, 1, 1],  # point 7
    ],
    np.double,
)

# these are the "midside" points of a quad cell.  See the definition of a
# vtkQuadraticHexahedron at:
# https://vtk.org/doc/nightly/html/classvtkQuadraticHexahedron.html
quad_pts = np.array(
    [
        (lin_pts[1] + lin_pts[0]) / 2,  # between point 0 and 1
        (lin_pts[1] + lin_pts[2]) / 2,  # between point 1 and 2
        (lin_pts[2] + lin_pts[3]) / 2,  # and so on...
        (lin_pts[3] + lin_pts[0]) / 2,
        (lin_pts[4] + lin_pts[5]) / 2,
        (lin_pts[5] + lin_pts[6]) / 2,
        (lin_pts[6] + lin_pts[7]) / 2,
        (lin_pts[7] + lin_pts[4]) / 2,
        (lin_pts[0] + lin_pts[4]) / 2,
        (lin_pts[1] + lin_pts[5]) / 2,
        (lin_pts[2] + lin_pts[6]) / 2,
        (lin_pts[3] + lin_pts[7]) / 2,
    ]
)

# introduce a minor variation to the location of the mid-side points
quad_pts += np.random.random(quad_pts.shape) * 0.3
pts = np.vstack((lin_pts, quad_pts))

# create the grid
cells = np.hstack((20, np.arange(20))).astype(np.int64, copy=False)
celltypes = np.array([CellType.QUADRATIC_HEXAHEDRON])
grid = pv.UnstructuredGrid(cells, celltypes, pts)
grid.plot(show_edges=True)
# finally, extract the surface and plot it
surf = grid.extract_surface()
surf.plot(show_scalar_bar=False)