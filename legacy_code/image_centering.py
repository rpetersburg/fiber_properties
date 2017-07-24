#==== Image Centering ====================================================#
#=========================================================================#

def set_fiber_data(self, method, **kwargs):
    """Set the fiber center, diameter, and centroid using the same method

    Args
    ----
    method : {'edge', 'radius', 'gaussian', 'circle'}
        Uses the respective method to find the fiber center
    **kwargs
        The keyworded arguments to pass to the centering method

    Sets
    ----
    _centroid.method : Pixel
        The centroid of the image in the context of the given method
    _center.method : Pixel
        The center of the fiber face in the context of the given method
    _diameter.method : float
        The diameter of the fiber face in the context of the given method
    """
    self.set_fiber_center(method, **kwargs)
    self.set_fiber_centroid(method, **kwargs)

def set_fiber_diameter(self, method, **kwargs):
    """Set the fiber diameter using given method

    Args
    ----
    method : {'edge', 'radius', 'gaussian', 'circle'}
        Uses the respective method to find the fiber center
    **kwargs :
        The keyworded arguments to pass to the centering method

    Sets
    ----
    _diameter.method : float
        The diameter of the fiber face in the context of the given method
    _center.method : Pixel
        The center of the fiber face in the context of the given method

    Raises
    ------
    RuntimeError
        cannot accept the 'circle' method when setting the diameter since
        it requires a known radius to run
    """
    if method == 'circle':
        raise RuntimeError('Fiber diameter cannot be set by circle method')
    self.set_fiber_center(method, **kwargs)

def set_fiber_center(self, method, show_image=False, **kwargs):
    """Find fiber center using given method

    Args
    ----
    method : {'edge', 'radius', 'gaussian', 'circle'}
        Uses the respective method to find the fiber center
    show_image : boolean, optional (default=False)
        Whether or not to show relevant fitting image
    **kwargs :
        The keyworded arguments to pass to the centering method

    Raises
    ------
    RuntimeError
        needs a valid method string to run the proper algorithm
    """
    center, diameter = fiber_center_and_diameter(self, method, show_image, **kwargs)
    setattr(self._center, method, center)
    setattr(self._diameter, method, diameter)

    # # Reset the fits due to new fiber parameters
    # if method == 'radius':
    #     self.set_fiber_center_radius_method(**kwargs)
    # elif method == 'edge':
    #     self.set_fiber_center_edge_method()
    # elif method == 'circle':
    #     self.set_fiber_center_circle_method(**kwargs)
    # elif method == 'gaussian':
    #     self.set_fiber_center_gaussian_method()
    # else:
    #     raise RuntimeError('Incorrect string for fiber centering method')

    # if show_image:
    #     center = getattr(self._center, method)
    #     r = getattr(self._diameter, method) / 2.0
    #     image = self.get_filtered_image()

    #     if method == 'gaussian':
    #         plot_overlaid_cross_sections(image, self.get_gaussian_fit(),
    #                                      center)
    #         plot_dot(image, center)
    #         show_plots()
    #     else:
    #         plot_image(remove_circle(image, center, r, res=1))
    #         plot_overlaid_cross_sections(image, image.max() / 2.0
    #                                      * circle_array(self.get_mesh_grid(),
    #                                                     center.x, center.y, r, res=1),
    #                                      center)
    #         if method == 'edge':
    #             for corner in self._edges:
    #                 plot_dot(image, corner)
    #         show_plots()

def set_fiber_center_gaussian_method(self):
    """Set fiber center using a Gaussian Fit

    Uses Scipy.optimize.curve_fit method to fit fiber image to
    gaussian_array(). The radius found extends to 2-sigma of the gaussian
    therefore encompassing ~95% of the imaged light. Use previous methods
    of center-finding to approximate the location of the center

    Sets
    ----
    _diameter.gaussian : float
        Diameter of the fiber in the gaussian method context
    _center.gaussian : {'x': float, 'y': float}
        Center of the fiber in the gaussian method context
    _fit.gaussian : 2D numpy.ndarray
        Best gaussian fit for the fiber image
    """
    _, coeffs = self.get_gaussian_fit(full_output=True)

    self._center.gaussian.x = coeffs[0]
    self._center.gaussian.y = coeffs[1]
    self._diameter.gaussian = abs(coeffs[2]) * 2.0
    self._gaussian_amp = coeffs[3]
    self._gaussian_offset = coeffs[4]

def set_fiber_center_radius_method(self, radius_tol=.03, radius_range=None, **kwargs):
    """Set fiber center using dark circle with varying radius

    Uses a golden mean optimization method to find the optimal radius of the
    dark circle that covers the fiber image used in
    get_fiber_centerCircleMethod(). The optimization is for a parameter
    array_sum which is weighted by the area of the circle, meaning that a
    smaller circle is preferred over one that simply covers the entire image

    Args
    ----
    radius_tol : number (default=1)
        Minimum possible range of radius values before ending iteration
    radius_range: int (in pixels)
        Range of tested radii, i.e. max(radius) - min(radius). If None,
        uses full possible range

    Sets
    ----
    _diameter.radius : float
        Diameter of the fiber in the radius method context
    _center.radius : {'x': float, 'y': float}
        Center of the fiber in the radius method context
    _diameter.circle : float
        Also uses the circle method, therefore changes this value
    _center.circle : float
        Also uses the circle method, therefore chnages this value
    """
    image = self.get_filtered_image()

    # Initialize range of tested radii
    r = np.zeros(4).astype(float)

    if radius_range is not None:
        approx_radius = self.get_fiber_radius(method='edge')
        radius_range /= 2.0

        r[0] = approx_radius - radius_range
        if r[0] < 0.0:
            r[0] = 0.0
        r[3] = approx_radius + radius_range
    else:
        r[0] = 0
        r[3] = min(self.height, self.width) / 2.0

    r[1] = r[0] + (1 - self._phi) * (r[3] - r[0])
    r[2] = r[0] + self._phi * (r[3] - r[0])

    array_sum = np.zeros(2).astype(float)
    for i in xrange(2):
        self.set_fiber_center(method='circle', radius=r[i+1],
                              image=image, **kwargs)
        array_sum[i] = (self._array_sum.circle
                        + self.threshold
                        * np.pi * r[i+1]**2)

    min_index = np.argmin(array_sum) # Integer 0 or 1 for min of r[1], r[2]

    while abs(r[3]-r[0]) > radius_tol:
        if min_index == 0:
            r[3] = r[2]
            r[2] = r[1]
            r[1] = r[0] + (1 - self._phi) * (r[3] - r[0])
        else:
            r[0] = r[1]
            r[1] = r[2]
            r[2] = r[0] + self._phi * (r[3] - r[0])

        array_sum[1 - min_index] = array_sum[min_index]

        self.set_fiber_center(method='circle', radius=r[min_index+1],
                              image=image, **kwargs)
        array_sum[min_index] = (self._array_sum.circle
                                + self.threshold
                                * np.pi * r[min_index+1]**2)

        min_index = np.argmin(array_sum) # Integer 0 or 1 for min of r[1], r[2]

    self._diameter.radius = r[min_index+1] * 2
    self._center.radius.y = self._center.circle.y
    self._center.radius.x = self._center.circle.x
    self._array_sum.radius = np.amin(array_sum)

def set_fiber_center_circle_method(self, radius=None, center_tol=.03,
                                   center_range=None, image=None, **kwargs):
    """Finds fiber center using a dark circle of set radius

    Uses golden mean method to find the optimal center for a circle
    covering the fiber image. The optimization is for a parameter array_sum
    that simply sums over the entire fiber image array

    Args
    ----
    radius : float
        Radius to use when creating circle
    center_tol : number (default=1)
        Minimum possible range of center values before ending iteration
    center_range: int (in pixels)
        Range of tested centers, i.e. max(x0) - min(x0). If None,
        uses full possible range
    image : 2d numpy.ndarray, optional
        The image being analyzed. This is only useful for the radius_method.
        Probably not for use outside the class.

    Sets
    ----
    _diameter.circle : float
        Diameter of the fiber in the circle method context
    _center.circle : {'x': float, 'y': float}
        Center of the fiber in the circle method context
    _diameter.edge : float
        If center_range is not None, approximates the circle's center using
        the edge method
    _center.edge : float
        If center_range is not None, approximates the circle's center using
        the edge method
    """
    res = int(1.0/center_tol)
    if image is None:
        image = self.get_filtered_image()
    if radius is None:
        radius = self.get_fiber_radius(method='edge')

    # Create four "corners" to test center of the removed circle
    x = np.zeros(4).astype(float)
    y = np.zeros(4).astype(float)

    if center_range is not None:
        approx_center = self.get_fiber_center(method='edge')
        center_range = center_range / 2.0

        x[0] = approx_center.x - center_range
        if x[0] < radius:
            x[0] = radius
        x[3] = approx_center.x + center_range
        if x[3] > self.width - radius:
            x[3] = self.width - radius

        y[0] = approx_center.y - center_range
        if y[0] < radius:
            y[0] = radius
        y[3] = approx_center.y + center_range
        if y[3] > self.height - radius:
            y[3] = self.height - radius

    else:
        x[0] = radius
        x[3] = self.width - radius

        y[0] = radius
        y[3] = self.height - radius

    x[1] = x[0] + (1 - self._phi) * (x[3] - x[0])
    x[2] = x[0] + self._phi * (x[3] - x[0])

    y[1] = y[0] + (1 - self._phi) * (y[3] - y[0])
    y[2] = y[0] + self._phi * (y[3] - y[0])

    # Initialize array sums to each corner
    array_sum = np.zeros((2, 2)).astype(float)
    for i in xrange(2):
        for j in xrange(2):
            removed_circle_array = remove_circle(image,
                                                 Pixel(x[i+1], y[j+1]),
                                                 radius, res=1)
            array_sum[j, i] = sum_array(removed_circle_array)

    # Find the index of the corner with minimum array_sum
    min_index = np.unravel_index(np.argmin(array_sum), (2, 2)) # Tuple

    while abs(x[3] - x[0]) > center_tol and abs(y[3] - y[0]) > center_tol:
        # Move the other corners to smaller search area
        if min_index[0] == 0:
            y[3] = y[2]
            y[2] = y[1]
            y[1] = y[0] + (1 - self._phi) * (y[3] - y[0])
        else:
            y[0] = y[1]
            y[1] = y[2]
            y[2] = y[0] + self._phi * (y[3] - y[0])
        if min_index[1] == 0:
            x[3] = x[2]
            x[2] = x[1]
            x[1] = x[0] + (1 - self._phi) * (x[3] - x[0])
        else:
            x[0] = x[1]
            x[1] = x[2]
            x[2] = x[0] + self._phi * (x[3] - x[0])

        # Replace the opposite corner array sum (so it doesn't need to be recalculated)
        array_sum[1 - min_index[0], 1 - min_index[1]] = array_sum[min_index]
        min_index = (1 - min_index[0], 1 - min_index[1])

        # Recalculate new sums for all four corners
        for i in xrange(2):
            for j in xrange(2):
                if i != min_index[1] or j != min_index[0]:
                    temp_res = 1
                    if abs(x[3] - x[0]) < 10*center_tol and abs(y[3] - y[0]) < 10*center_tol:
                        temp_res = res
                    removed_circle_array = remove_circle(image,
                                                         Pixel(x[i+1], y[j+1]),
                                                         radius, temp_res)
                    array_sum[j, i] = sum_array(removed_circle_array)

        min_index = np.unravel_index(np.argmin(array_sum), (2, 2))

    self._center.circle.x = x[min_index[1]+1]
    self._center.circle.y = y[min_index[0]+1]
    self._diameter.circle = radius * 2.0
    self._array_sum.circle = np.amin(array_sum)

def set_fiber_center_edge_method(self, **kwargs):
    """TAverages the fiber edges to set the fiber center

    Sets
    ----
    self._center.edge.y : float
    self._center.edge.x : float
    """
    self.set_fiber_edges(**kwargs)

    self._center.edge.y = (self._edges.top.y + self._edges.bottom.y) / 2.0
    self._center.edge.x = (self._edges.left.x + self._edges.right.x) / 2.0