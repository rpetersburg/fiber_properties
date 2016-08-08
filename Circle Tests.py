from NumpyArrayHandler import NumpyArrayHandler as NAH
import numpy as np
import time

length = 1000

mesh_grid = np.meshgrid(np.arange(length), np.arange(length))

errors = []
times = [] 
for i in xrange(30):
    t0 = time.time()

    radius = (np.random.rand()*0.1 + 0.2)*length
    x0 = np.random.rand()+length/2.0
    y0 = np.random.rand()+length/2.0
    #radius = length * (0.2 + 0.01 * i)
    #x0 = length * 0.5 + 0.43
    #y0 = length * 0.5 + 0.76

    area = np.pi * radius**2

    circle = NAH.circleArray(mesh_grid, x0, y0, radius, res=1)
    #NAH.showImageArray(circle)
    
    errors.append((circle.sum() - area)/area*100)
    #print 'x0:', x, 'y0:', y

    t1 = time.time()
    times.append(t1-t0)

    print 'Error:', errors[i], '%  Radius:', radius, '  Time:', t1-t0

errors = np.array(errors)
print "Error mean:", np.abs(errors).mean(), '%'
print "Error var:", errors.var()

times = np.array(times)
print "Time mean:", times.mean()