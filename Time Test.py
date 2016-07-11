import time
import numpy as np

a, b = np.meshgrid(np.linspace(0,1000,1001),
                   np.linspace(0,1000,1001))

res = 1000
x = 1
y = 1
vals = np.arange(0,1,1.0/res) - 0.5 + 0.5/res

t0 = time.time()

for x_temp in vals:
    for y_temp in vals:
        xy = x_temp**2 + y_temp**2

t1 = time.time()

t2 = time.time()

for i in xrange(res):
    for j in xrange(res):
        x_temp = x - 0.5 + (i+0.5) / float(res)
        y_temp = y - 0.5 + (j+0.5) / float(res)
        xy = x_temp**2 + y_temp**2

t3 = time.time()

print t1-t0
print t3-t2