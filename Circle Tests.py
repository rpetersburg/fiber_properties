import NumpyArrayHandler as NAH
import numpy as np

obj = NAH.NumpyArrayHandler()

mesh_grid = np.meshgrid(np.linspace(0, 1000, 1000),
                        np.linspace(0, 1000, 1000))

errors = []  
for i in xrange(100):
    radius = np.random.rand()*300
    x = np.random.rand()+500
    y = np.random.rand()+500
    #print radius, x, y
    circle = obj.circleArray(mesh_grid, x, y, radius)
    area = np.pi * radius**2
    print 'Error:', np.abs(np.sum(circle) - area)
    errors.append(np.abs(np.sum(circle) - area))
errors = np.array(errors)
print "Error mean:", errors.mean()