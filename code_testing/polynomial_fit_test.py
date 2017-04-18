from fiber_properties import polynomial_fit, polynomial_array, FiberImage
import numpy as np
import matplotlib.pyplot as plt

image = np.meshgrid(np.arange(100), np.arange(100))
coeffs = np.random.rand(6)
print coeffs
image = polynomial_array(image, *coeffs).reshape(100,100).astype('float64')
# image += np.random.rand(100,100)*10.0 - 5.0
plt.figure(1)
plt.imshow(image)

fit = polynomial_fit(image, 2)
plt.figure(2)
plt.imshow(fit)

plt.show()