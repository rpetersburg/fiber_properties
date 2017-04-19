if self.camera == 'in':
    fiber_face = circle_array(self.get_mesh_grid(), center.x, center.y, radius)
    fiber_face *= np.median(intensity_array(image, center, radius))
    gaussian_fit += fiber_face




center = self.get_fiber_center(method='edge')
radius = self.get_fiber_radius(method='edge')
image = self.get_filtered_image()

if self.camera == 'in':
    test_radius = -1
    factor = 0.9
    fiber_face = circle_array(self.get_mesh_grid(), center.x,
                              center.y, radius, res=1)
    fiber_face *= np.median(intensity_array(image, center, radius))
    while test_radius < 1 or test_radius > radius:
        factor += 0.1
        initial_guess = (center.x, center.y, 100 / self.get_pixel_size(),
                         image.max(), image.min())
        try:
            fit, opt_parameters = gaussian_fit(image,
                                               initial_guess=initial_guess,
                                               full_output=True,
                                               center=center,
                                               radius=radius)
            test_radius = abs(opt_parameters[2])
        except RuntimeError:
            test_radius = -1

    x0 = opt_parameters[0]
    y0 = opt_parameters[1]
    amp = opt_parameters[3]
    offset = opt_parameters[4]

else:
    initial_guess = (fiber_center.x, fiber_center.y, fiber_radius,
                     filtered_image.max(), filtered_image.min())

    _, opt_parameters = gaussian_fit(filtered_image,
                                     initial_guess=initial_guess,
                                     full_output=True,
                                     center=center,
                                     radius=radius)
    x0 = opt_parameters[0]
    y0 = opt_parameters[1]
    test_radius = abs(opt_parameters[2])
    amp = opt_parameters[3]
    offset = opt_parameters[4]

self._center.gaussian.x = x0
self._center.gaussian.y = y0
self._diameter.gaussian = test_radius * 2.0
self._gaussian_amp = amp
self._gaussian_offset = offset