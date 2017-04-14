import time
from scipy.signal import medfilt2d, medfilt, order_filter
import numpy as np
import matplotlib.pyplot as plt

a = 2**16 * np.random.rand(100,100)
plt.figure(1)
plt.imshow(a)

start = time.time()
b = medfilt2d(a, 9)
end = time.time()
print(end-start)

plt.figure(2)
plt.imshow(b)

start = time.time()
c = medfilt2d(medfilt2d(a, 9), 9)
end = time.time()
print(end-start)

plt.figure(3)
plt.imshow(c)

start = time.time()
d = medfilt2d(a, 19)
end = time.time()
print(end-start)

plt.figure(4)
plt.imshow(d)

plt.show()
