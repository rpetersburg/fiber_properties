from fiber_properties import rectangle_array
import numpy as np
import matplotlib.pyplot as plt

a = np.meshgrid(np.arange(512), np.arange(512))
rect = rectangle_array(a, 200, 300, 440, 110, 10)

plt.imshow(rect)
plt.show()
