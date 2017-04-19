from fiber_properties import polynomial_fit, polynomial_array, FiberImage, crop_image, plot_overlaid_cross_sections
import numpy as np
import matplotlib.pyplot as plt

# image = np.meshgrid(np.arange(100), np.arange(100))
# coeffs = np.random.rand(10) * 10.0
# print coeffs
# image = polynomial_array(image, *coeffs).reshape(100,100).astype('float64')
# image += np.random.rand(100,100)*10.0 - 5.0

im_obj = FiberImage('../data/modal_noise/Kris_data/circular_100um/coupled_agitation/nf_corrected.fit')
image = im_obj.get_image()
center = im_obj.get_fiber_center()
# radius = im_obj.get_fiber_radius()
# image, center = crop_image(image, center=center, radius=radius/2)

plt.figure(1)
plt.imshow(image)

# fit = polynomial_fit(image, 6)
fit = im_obj.get_polynomial_fit()
plt.figure(2)
plt.imshow(fit)

plot_overlaid_cross_sections(image, fit, center)

plt.show()