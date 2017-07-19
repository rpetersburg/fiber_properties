"""fiber_center.py was written by Ryan Petersburg for use with fiber
characterization on the EXtreme PREcision Spectrograph

The functions in this module calculate the center and various dimensions of a
FiberImage object
"""
from .numpy_array_handler import circle_array, remove_circle, sum_array, mesh_grid_from_array
from .plotting import plot_image, plot_dot, plot_overlaid_cross_sections, show_plots
from .containers import Pixel
from .fiber_centroid import calc_centroid
import numpy as np

# Golden ratio for the optimization tests
PHI = (5 ** 0.5 - 1) / 2

def fiber_center_and_diameter(im_obj, method, show_image=False, **kwargs):
    """Find fiber center and diameter using the given method

    Args
    ----
    method : {'edge', 'radius', 'gaussian', 'circle'}
        Uses the respective method to find the fiber center
    show_image : boolean, optional (default=False)
        Whether or not to show relevant fitting images
    **kwargs :
        The keyworded arguments to pass to the centering method

    Raises
    ------
    RuntimeError
        needs a valid method string to run the proper algorithm
    """
    if method == 'radius':
        center, diameter = _radius_method(im_obj, **kwargs)
    elif method == 'edge':
        center, diameter = _edge_method(im_obj, **kwargs)
    elif method == 'circle':
        center, diameter = _circle_method(im_obj, **kwargs)
    elif method == 'gaussian':
        center, diameter = _gaussian_method(im_obj, **kwargs)
    elif method == 'full':
        center, diameter = _full_method(im_obj, **kwargs)
    else:
        raise RuntimeError('Incorrect string for fiber centering method')

    if show_image:
        radius = diameter / 2.0
        image = im_obj.get_image()

        if method == 'gaussian':
            plot_overlaid_cross_sections(image, 
                                         im_obj.get_gaussian_fit(),
                                         center)
            plot_dot(image, center)
        else:
            plot_image(remove_circle(image, center, radius, res=1))
            plot_overlaid_cross_sections(image,
                                         circle_array(im_obj.get_mesh_grid(),
                                                        center.x, center.y,
                                                        radius, res=1)
                                         *image.max() / 2.0,
                                         center)
            plot_dot(image, center)
            if method == 'edge':
                for corner in im_obj._edges:
                    plot_dot(image, corner)
        show_plots()

    return center, diameter

def _full_method(im_obj, kernel=None, threshold=None, **kwargs):
    """Centroids a boolean image above the FiberImage threshold.

    Returns
    -------
    center : Pixel
    diameter : float (pixels)
    """
    if threshold is None:
        threshold = im_obj.threshold
    image = (im_obj.get_filtered_image(kernel) > threshold).astype('uint8')
    center = calc_centroid(image)
    x_array, y_array = mesh_grid_from_array(image)
    dist_image = np.sqrt((x_array - center.x)**2 + (y_array - center.y)**2)
    dist_image *= image
    diameter = dist_image.max() * 2.0
    # _, diameter = _edge_method(im_obj, kernel=kernel, **kwargs)
    return center, diameter

def _edge_method(im_obj, **kwargs):
    """Averages the fiber edges to set the fiber center

    Returns
    -------
    center : Pixel
    diameter : float (pixels)
    """
    im_obj.set_fiber_edges(**kwargs)
    edges = im_obj._edges

    y = (edges.top.y + edges.bottom.y) / 2.0
    x = (edges.left.x + edges.right.x) / 2.0
    center = Pixel(x,y)

    # average the horizontal and vertical distances
    diameter = (np.sqrt(((edges.right - edges.left)**2).as_array().sum())
                + np.sqrt(((edges.bottom - edges.top)**2).as_array().sum())) / 2.0

    return center, diameter

def _gaussian_method(im_obj, **kwargs):
    """Set fiber center using a Gaussian Fit

    Uses Scipy.optimize.curve_fit method to fit fiber image to
    gaussian_array(). The radius found extends to 2-sigma of the gaussian
    therefore encompassing ~95% of the imaged light. Use previous methods
    of center-finding to approximate the location of the center

    Returns
    -------
    center : Pixel
        Center of the fiber in the gaussian method context
    diameter : float (pixels)
    """
    if im_obj.gaussian_coeffs is None:
        im_obj.set_gaussian_fit(**kwargs)
    coeffs = im_obj.gaussian_coeffs
    center = Pixel(coeffs[0], coeffs[1])
    diameter = abs(coeffs[2]) * 2.0
    return center, diameter

def _radius_method(im_obj, radius_tol=.03, radius_range=None, option=None,
                   kernel=None, threshold=None, **kwargs):
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

    Returns
    -------
    center : Pixel
    diameter : float (pixels)
    array_sum : float
        If option is 'all'
    """
    image = im_obj.get_filtered_image(kernel)
    if threshold is None:
        threshold = im_obj.threshold

    # Initialize range of tested radii
    r = np.zeros(4).astype(float)

    if radius_range is not None:
        approx_radius = im_obj.get_fiber_radius(method='edge')
        radius_range /= 2.0

        r[0] = approx_radius - radius_range
        if r[0] < 0.0:
            r[0] = 0.0
        r[3] = approx_radius + radius_range
    else:
        r[0] = 0
        r[3] = min(im_obj.height, im_obj.width) / 2.0

    r[1] = r[0] + (1 - PHI) * (r[3] - r[0])
    r[2] = r[0] + PHI * (r[3] - r[0])

    array_sum = np.zeros(2).astype(float)
    for i in xrange(2):
        center, _, array_sum[i] = _circle_method(im_obj, image=image,
                                                 radius=r[i+1],
                                                 option='all', **kwargs)
        array_sum[i] += im_obj.threshold * np.pi * r[i+1]**2

    min_index = np.argmin(array_sum) # Integer 0 or 1 for min of r[1], r[2]

    while abs(r[3]-r[0]) > radius_tol:
        if min_index == 0:
            r[3] = r[2]
            r[2] = r[1]
            r[1] = r[0] + (1 - PHI) * (r[3] - r[0])
        else:
            r[0] = r[1]
            r[1] = r[2]
            r[2] = r[0] + PHI * (r[3] - r[0])

        array_sum[1 - min_index] = array_sum[min_index]

        center, _, array_sum[min_index] = _circle_method(im_obj, image=image,
                                                         radius=r[min_index+1],
                                                         option='all', **kwargs)
        array_sum[min_index] += threshold * np.pi * r[min_index+1]**2

        min_index = np.argmin(array_sum) # Integer 0 or 1 for min of r[1], r[2]

    array_sum = np.amin(array_sum)
    diameter = r[min_index+1] * 2

    if option == 'all':
        return center, diameter, array_sum
    return center, diameter

def _circle_method(im_obj, image=None, radius=None, center_tol=.03,
                  center_range=None, option=None, kernel=None, **kwargs):
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

    Returns
    -------
    center : Pixel
    diameter : float (pixels)
    array_sum : float
        if option is 'all'
    """
    res = int(1.0/center_tol)
    if image is None:
        image = im_obj.get_filtered_image(kernel)
    if radius is None:
        radius = im_obj.get_fiber_radius(method='edge')

    # Create four "corners" to test center of the removed circle
    x = np.zeros(4).astype(float)
    y = np.zeros(4).astype(float)

    if center_range is not None:
        approx_center = im_obj.get_fiber_center(method='edge')
        center_range = center_range / 2.0

        x[0] = approx_center.x - center_range
        if x[0] < radius:
            x[0] = radius
        x[3] = approx_center.x + center_range
        if x[3] > im_obj.width - radius:
            x[3] = im_obj.width - radius

        y[0] = approx_center.y - center_range
        if y[0] < radius:
            y[0] = radius
        y[3] = approx_center.y + center_range
        if y[3] > im_obj.height - radius:
            y[3] = im_obj.height - radius

    else:
        x[0] = radius
        x[3] = im_obj.width - radius

        y[0] = radius
        y[3] = im_obj.height - radius

    x[1] = x[0] + (1 - PHI) * (x[3] - x[0])
    x[2] = x[0] + PHI * (x[3] - x[0])

    y[1] = y[0] + (1 - PHI) * (y[3] - y[0])
    y[2] = y[0] + PHI * (y[3] - y[0])

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
            y[1] = y[0] + (1 - PHI) * (y[3] - y[0])
        else:
            y[0] = y[1]
            y[1] = y[2]
            y[2] = y[0] + PHI * (y[3] - y[0])
        if min_index[1] == 0:
            x[3] = x[2]
            x[2] = x[1]
            x[1] = x[0] + (1 - PHI) * (x[3] - x[0])
        else:
            x[0] = x[1]
            x[1] = x[2]
            x[2] = x[0] + PHI * (x[3] - x[0])

        # Replace the opposite corner array sum (so it doesn't need to be recalculated)
        array_sum[1 - min_index[0], 1 - min_index[1]] = array_sum[min_index]
        min_index = (1 - min_index[0], 1 - min_index[1])

        # Recalculate new sums for all four corners
        for i in xrange(2):
            for j in xrange(2):
                if i != min_index[1] or j != min_index[0]:
                    temp_res = 1
                    if (abs(x[3] - x[0]) < 10*center_tol
                            and abs(y[3] - y[0]) < 10*center_tol):
                        temp_res = res
                    removed_circle_array = remove_circle(image,
                                                         Pixel(x[i+1], y[j+1]),
                                                         radius, temp_res)
                    array_sum[j, i] = sum_array(removed_circle_array)

        min_index = np.unravel_index(np.argmin(array_sum), (2, 2))

    center = Pixel(x[min_index[1]+1], y[min_index[0]+1])
    diameter = radius * 2
    array_sum = np.amin(array_sum)

    if option == 'all':
        return center, diameter, array_sum
    return center, diameter
