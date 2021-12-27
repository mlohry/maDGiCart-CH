import pyvista as pv
import sys

# read the data
filename = sys.argv[1]
grid = pv.read(filename)

# plot the data with an automatically created Plotter
grid.contour()

print(f"All arrays: {grid.array_names}")
plotter = pv.Plotter()
plotter.add_mesh(grid, scalars='c', cmap='inferno')
plotter.view_xy()
plotter.show()
