# %%
import pyvista as pv
from pyvistaqt import BackgroundPlotter


# %%
sphere = pv.Sphere()
# %%

plotter = BackgroundPlotter()
# %%
plotter.add_mesh(sphere)
# %%
plotter.add_mesh(sphere, color='r',style='wireframe')
# %%
waiting = input("Press Enter to continue.")

