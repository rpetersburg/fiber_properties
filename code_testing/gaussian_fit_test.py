from fiber_properties import gaussian_fit, gaussian_array, FiberImage
import numpy as np
import matplotlib.pyplot as plt

# image = np.meshgrid(np.arange(100), np.arange(100))
# coeffs = np.random.normal(size=5) * 15 + 50
# print coeffs
# image = gaussian_array(image, *coeffs).reshape(100,100).astype('float64')
# image += np.random.normal(size=image.shape)

im_obj = FiberImage('../data/modal_noise/Kris_data/rectangular_100x300um/coupled_agitation/ff_corrected.fit')
image = im_obj.get_image()
center = im_obj.get_fiber_center()

plt.figure(1)
plt.imshow(image)

fit, output = gaussian_fit(image, full_output=True)
print output
plt.figure(2)
plt.imshow(fit)

plt.show()
