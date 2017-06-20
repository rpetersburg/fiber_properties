from fiber_properties import FiberImage, circle_array
import numpy as np

mesh_grid = np.meshgrid(np.arange(1000), np.arange(1000))
base = 50000.0 * circle_array(mesh_grid, 500, 500, 50, res=16) + 400.0
dark = np.random.normal(base, 5.0*np.sqrt(base)).astype('float64')
# dark = np.random.poisson(25.0*base).astype('float64') / 25.0
FiberImage(dark).save_image('noise_test.fit')