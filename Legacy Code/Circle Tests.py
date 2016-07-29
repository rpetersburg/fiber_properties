import NumpyArrayHandler as NAH
import numpy as np
import time

obj = NAH.NumpyArrayHandler()

length = 1000

mesh_grid = np.meshgrid(np.linspace(0, length, length+1),
                        np.linspace(0, length, length+1))

errors = []
times = [] 
for i in xrange(50):
    t0 = time.time()

    radius = (np.random.rand()*0.1 + 0.2)*length
    x = np.random.rand()+length/2.0
    y = np.random.rand()+length/2.0
    #radius = 200 + 10 * i
    #x = 500.43
    #y = 500.76

    area = np.pi * radius**2

    circle = obj.circleArray(mesh_grid, x, y, radius)
    
    errors.append((np.sum(circle) - area)/area*100)
    #print 'x0:', x, 'y0:', y
    print 'Error:', errors[i], '%     Radius:', radius

    t1 = time.time()
    times.append(t1-t0)

errors = np.array(errors)
print "Error mean:", np.abs(errors).mean(), '%'
print "Error var:", errors.var()

times = np.array(times)
print "Time mean:", times.mean()