from fiber_properties.numpy_array_handler import _median
import numpy as np

a = np.random.rand(11)
print _median(a)
print np.sort(a)[int((a.size-1)/2)]
print np.sort(a)