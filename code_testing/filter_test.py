import time
from scipy.signal import medfilt2d, medfilt, order_filter
import numpy as np

a = 2**16 * np.random.rand(100,100)

start = time.time()
b = medfilt2d(a, 9)
end = time.time()
print(end-start)

a = 2**16 * np.random.rand(200,200)

start = time.time()
b = medfilt2d(a, 9)
end = time.time()
print(end-start)

a = 2**16 * np.random.rand(400,400)

start = time.time()
b = medfilt2d(a, 9)
end = time.time()
print(end-start)

a = 2**16 * np.random.rand(800,800)

start = time.time()
b = medfilt2d(a, 9)
end = time.time()
print(end-start)
