import time
import numpy as np

t0 = time.time()

for i in range(10000):
    for j in range(10000):
        pass

t1 = time.time()

t2 = time.time()

for i in range(10000):
    for j in range(10000):
        a = i**2 + j**2

t3 = time.time()

print t1-t0
print t3-t2