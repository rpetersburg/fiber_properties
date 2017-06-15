"""NumpyArrayHandler.py was written by Ryan Petersburg for use with fiber
characterization on the EXtreme PREcision Spectrograph

The functions in this module are used to handle two dimensional np.ndarray
objects that represent images from optical fibers. Functions include
array summing, image cropping, image filtering, function fitting, and
function generation.
"""
import numpy as np
from scipy import optimize as opt
from scipy.signal import medfilt2d, order_filter
from PIL import Image, ImageDraw
from containers import Pixel
import math
from .filter_image import median, c_filter_image, c_filter_image_zero_fill

#=============================================================================#
#===== Array Summing =========================================================#
#=============================================================================#

def sum_array(image):
    """Returns the sum of all elements in a numpy.ndarray"""
    return np.sum(image)

def sum_rows(image):
    """Sums the rows of a 2D np.ndarray

    Args
    ----
    image : 2D numpy.ndarray

    Returns
    -------
    summed_row : 1D numpy.ndarray

    """
    row_sum = np.sum(image, axis=0)
    return ((row_sum - np.min(row_sum)) / image.shape[0]).astype('float64')

def sum_columns(image):
    """Sums the columnns of a 2D np.ndarray

    Args
    ----
    image : 2D numpy.ndarray

    Returns
    -------
    summed_column : 1D numpy.ndarray
    """
    column_sum = np.sum(image, axis=1)
    return ((column_sum - np.min(column_sum)) / image.shape[1]).astype('float64')

#=============================================================================#
#===== Array Alterations =====================================================#
#=============================================================================#

def mesh_grid_from_array(image):
    """Creates a numpy meshgrid of pixel number for an image

    Args
    ----
    image : 2D numpy.ndarray

    Returns
    -------
    mesh_grid_x : 2D numpy.ndarray
        The x position of each point in the grid
    mesh_grid_y : 2D numpy.ndarray
        The y position of each point in the grid
    """
    return np.meshgrid(np.arange(image.shape[1]).astype('float64'),
                       np.arange(image.shape[0]).astype('float64'))

def intensity_array(image, center, radius):
    """Returns intensities from inside a circle

    Returns an array of intensities from image which are contained
    within the circle with radius centered at (x0, y0)

    Args
    ----
    image : 2D numpy.ndarray
    center : Pixel
    radius : number (pixels)

    Returns
    -------
    intensity_array : 1D numpy.ndarray
        Intensities of the elements contained within the given circle
    """
    image_crop, new_center = crop_image(image, center, radius)
    x, y = mesh_grid_from_array(image_crop)
    mask = (new_center.x-x)**2 + (new_center.y-y)**2 <= (radius)**2
    intensity_array = image_crop[mask]

    # height, width = image_crop.shape
    # intensity_list = []
    # for x in xrange(width):
    #     for y in xrange(height):
    #         if (new_center.x-x)**2 + (new_center.y-y)**2 <= (radius)**2:
    #             intensity_list.append(image_crop[y, x])
    # intensity_array = np.array(intensity_list)

    return intensity_array

def crop_image(image, center, radius, full_output=True):
    """Crops image to square with radius centered at (y0, x0)

    Args
    ----
    image : 2D numpy.ndarray
    center : Pixel
    radius : number (pixels)

    Returns
    -------
    image_crop : 2D numpy.ndarray
    new_center : Pixel
    """
    top = int(center.y-radius)
    if top < 0:
        top = 0
    left = int(center.x-radius)
    if left < 0:
        left = 0
    if (int(center.y) == center.y and int(center.x) == center.x
                                  and int(radius) == radius):
        image_crop = image[top : int(center.y + radius) + 1,
                           left : int(center.x + radius) + 1]
    else:
        image_crop = image[top : int(center.y + radius) + 2,
                           left : int(center.x + radius) + 2]

    if full_output:
        new_center = Pixel(center.x - left, center.y - top)
        return image_crop, new_center
    return image_crop

def subframe_image(image, subframe_x, subframe_y, width, height):
    """Creates the subframe of an image with the given parameters."""
    return image[subframe_y : subframe_y + height,
                 subframe_x : subframe_x + width]

def remove_circle(image, center, radius, res=1):
    """Removes a circle from an array

    Args
    ----
    image : 2D numpy.ndarray
    center : Pixel
    radius : number (pixels)

    Returns
    -------
    removed_circle_array : 2D numpy.ndarray
        Input image array with the defined circle removed
    """
    mesh_grid = mesh_grid_from_array(image)
    return image * (1 - circle_array(mesh_grid, center.x,
                                     center.y, radius, res))

def isolate_circle(image, center, radius, res=1):
    """Isolates a circle in an array

    Args
    ----
    image : 2D numpy.ndarray
    center : Pixel
    radius : number (pixels)

    Returns
    -------
    isolated_circle_array : 2D numpy.ndarray
        Input image array with the defined circle isolated in the image
    """
    mesh_grid = mesh_grid_from_array(image)
    return image * circle_array(mesh_grid, center.x, center.y, radius, res)

def isolate_rectangle(image, corners=None, center=None, **kwargs):
    """Isolates a rectangle in an array

    Args
    ----
    image : 2D numpy.ndarray
    corners : list of Pixels, optional
    center : Pixel, optional

    Returns
    -------
    isolated_rectangle_array : 2D numpy.ndarray
        Input image array with the defined rectangle isolated in the image
    """
    mesh_grid = mesh_grid_from_array(image)
    if corners is not None:
        return image * rectangle_array(mesh_grid, corners=corners)
    return image * rectangle_array(mesh_grid, center.x, center.y, **kwargs)

def apply_window(image):
    """Applies a FFT window to an image

    Args
    ----
    image : 2D numpy.ndarray

    Returns
    -------
    windowed_array : 2D numpy.ndarray
    """
    height, width = image.shape
    x_array, y_array = mesh_grid_from_array(image)
    x0 = width/2
    y0 = height/2
    r_array = np.sqrt((x_array-x0)**2 + (y_array-y0)**2) + min(height, width) / 2
    window = hann_poisson_window(min(height, width), r_array)
    return image * window

def hann_poisson_window(arr_len, arr=None):
    """Hann-Poisson FFT window:
    https://en.wikipedia.org/wiki/Window_function#Hann.E2.80.93Poisson_window

    Args
    ----
    arr_len : int
        Length of the array on which the window is built
    arr : numpy.ndarray, optional (default=np.arange(arr_len))
        Any dimensional numpy array on which the window is applied

    Returns
    -------
    hann_poisson_window : 1D numpy.ndarray
    """
    if arr is None:
        arr = np.arange(arr_len)
    hann = 0.5 * (1 - np.cos(2*np.pi*arr / (arr_len - 1)))
    alpha = 2.0
    poisson = np.exp(-alpha/(arr_len-1) * np.abs(arr_len - 1 - 2*arr))
    return hann * poisson

def poisson_window(arr_len, arr=None):
    """Poisson FFT window:
    https://en.wikipedia.org/wiki/Window_function#Exponential_or_Poisson_window

    Args
    ----
    arr_len : int
        Length of the array on which the window is built
    arr : numpy.ndarray, optional (default=np.arange(arr_len))
        Any dimensional numpy array on which the window is applied

    Returns
    -------
    poisson_window : 1D numpy.ndarray

    """
    if arr is None:
        arr = np.arange(arr_len)
    tau = (arr_len / 2) * (8.69 / 60)
    poisson = np.exp(-np.abs(arr - (arr_len-1)/2) / tau)
    return poisson

def filter_image(image, kernel_size, quick=None, cython=False, zero_fill=False):
    """
    Args
    ----
    image : 2D numpy.ndarray
    kernel_size : int (odd)
        side length of the kernel which the median filter uses

    Returns
    -------
    filtered_image : 2D numpy.ndarray

    """
    if kernel_size < 2.0:
        return image
    height, width = image.shape

    if kernel_size > width or kernel_size > height:
        raise ValueError('Please choose kernel size smaller than image '
                         + 'dimensions')
    if kernel_size % 2.0 != 1.0:
        raise ValueError('Please use odd integer for kernel size')
    if kernel_size > 101:
        kernel_size = 101

    if quick is None:
        quick = False
        if kernel_size < 11:
            quick = True

    if quick:
        return medfilt2d(image, kernel_size)
    if cython:
        if zero_fill:
            return c_filter_image_zero_fill(image, kernel_size)
        return c_filter_image(image, kernel_size)

    radius = (kernel_size-1) / 2
    x_array, y_array = np.meshgrid(np.arange(kernel_size),
                                   np.arange(kernel_size))
    mask = (x_array-radius)**2 + (y_array-radius)**2 <= radius**2

    filtered_image = np.zeros_like(image)
    if zero_fill:
        zero_image = np.zeros([kernel_size, kernel_size])
        temp_mask = mask
    for y in xrange(height):
        y_edge = False
        if y < radius:
            mask_top = kernel_size - radius - 1 - y
            top = 0
            y_edge = True
        else:
            mask_top = 0
            top = y - radius
        if y > height - 1 - radius:
            mask_bottom = kernel_size - radius + (height - 1 - y)
            bottom = height
            y_edge = True
        else:
            mask_bottom = kernel_size
            bottom = y + radius + 1

        for x in xrange(width):
            x_edge = False
            if x < radius:
                mask_left = kernel_size - radius - 1 - x
                left = 0
                x_edge = True
            else:
                mask_left = 0
                left = x - radius
            if x > width - 1 - radius:
                mask_right = kernel_size - radius + (width - 1 - x)
                right = width
                x_edge = True
            else:
                mask_right = kernel_size
                right = x + radius + 1

            if zero_fill and (y_edge or x_edge):
                image_crop = zero_image.copy()
                image_crop[mask_top:mask_bottom, mask_left:mask_right] = image[top:bottom, left:right]
            else:
                image_crop = image[top:bottom, left:right]
            if not zero_fill:                
                temp_mask = mask[mask_top:mask_bottom, mask_left:mask_right]
            filtered_image[y, x] = median(image_crop[temp_mask])

    return filtered_image

#=============================================================================#
#===== 2D Array Functions ====================================================#
#=============================================================================#

def gaussian_array(mesh_grid, x0, y0, radius, amp, offset):
    """Creates a 2D gaussian function as a 1D array

    Args
    ----
    mesh_grid : numpy.meshgrid
    image : 2D numpy.ndarray
    x0 : number (pixels)
    y0 : number (pixels)
    radius : number (pixels)
        "Radius" of gaussian (2 standard deviations or 95% integrated)
    amp: number
        Amplitude of the gaussian

    Returns
    -------
    gaussian_array : 1D numpy.ndarray
        Ravelled gaussian numpy array usable in Scipy.optimize.curve_fit method
        and properly reshaped by gaussian_array.reshape(height, width)

    """
    gaussian_array = offset + amp * np.exp(-2*(mesh_grid[0] - x0)**2 / radius**2
                                           -2*(mesh_grid[1] - y0)**2 / radius**2)
    return gaussian_array.ravel()

def circle_array(mesh_grid, x0, y0, radius, res=1):
    """Creates a 2D tophat function of amplitude 1.0

    Args
    ----
    mesh_grid : numpy.meshgrid
    x0 : number (pixels)
    y0 : number (pixels)
    radius : number (pixels)
    res : int, optional (default=1)
        Resolution element for more precise circular edges

    Returns
    -------
    circle_array : 2D numpy.ndarray
        Points inside the circle are 1.0 and outside the circle are 0.0. Points
        along the edge of the circle are weighted based on their relative
        distance to the center

    """
    x0 = float(x0)
    y0 = float(y0)
    radius = float(radius)

    x_array = mesh_grid[0].astype('float64')
    y_array = mesh_grid[1].astype('float64')

    if float(res) <= 1.0:
        circle_array = ((x_array-x0)**2 + (y_array-y0)**2 <= radius**2).astype('float64')

    else:
        circle_array = ((x_array-x0)**2 + (y_array-y0)**2
                        < (radius - np.sqrt(2) / 2.0)**2).astype('float64')

        res_array = np.arange(-0.5, 0.5, 1.0 / res) + 0.5 / res
        res_mesh_x, res_mesh_y = np.meshgrid(res_array, res_array)
        res_val = 1.0 / res**2
        height, width = circle_array.shape

        for x in range(max(0, int(x0-radius)), min(width, int(x0+radius) + 2)):
            for y in range(max(0, int(y0-radius)), min(height, int(y0+radius) + 2)):
                if circle_array[y, x] < 1.0:
                    if (x-x0)**2 + (y-y0)**2 <= (radius + np.sqrt(2) / 2.0)**2:
                        circle_array[y, x] += res_val * ((res_mesh_x+x-x0)**2
                                                         + (res_mesh_y+y-y0)**2
                                                         <= radius**2
                                                        ).astype('float64').sum()
    return circle_array

def rectangle_array(mesh_grid, x0=None, y0=None, width=None, height=None, angle=None, corners=None):
    """Creates a 2D rectangle array of amplitude 1.0

    Args
    ----
    mesh_grid: numpy.meshgrid
    x0 : number (pixels)
    y0 : number (pixels)
    height : number (pixels)
    width : number (pixels)
    angle : number (degrees)
    corner : sequence of Pixel objects (pixels)

    Returns
    -------
    rectangle_array : 2D numpy.ndarray
        Points inside the rectangle are 1.0 and outside the rectangle are 0.0
    """
    image = Image.fromarray(np.zeros_like(mesh_grid[0]).astype('float64'))
    if corners is not None:
        ImageDraw.Draw(image).polygon([corner.as_tuple() for corner in corners], fill=1.0)
    else:
        x0 = float(x0)
        y0 = float(y0)
        width = float(width)
        height = float(height)
        angle = float(angle)

        rect = np.array([(0,0), (width, 0), (width, height), (0, height), (0,0)])
        theta = (np.pi / 180.0) * angle
        R = np.array([[np.cos(theta), -np.sin(theta)],
                      [np.sin(theta), np.cos(theta)]])
        rect = np.dot(rect, R)
        offset = np.array([x0-rect[2,0]/2.0, y0-rect[2,1]/2.0])
        rect += offset

        ImageDraw.Draw(image).polygon([tuple(p) for p in rect], fill=1.0)

    rectangle_array = np.asarray(image)
    return rectangle_array

def polynomial_array(mesh_grid, *coeffs):
    """2D polynomial of arbitrary degree for given x, y

    Uses a mesh grid and list of coefficients to create a two dimensional
    polynomial array in the following order:
    c0 + c1*x + c2*y + c3*x^2 + c4*x*y + c5*y^2 + c6*x^3 + c7*x^2*y + ...

    Args
    ----
    *coeffs :
        coefficients where the length is determined by the degree of the
        polynomial (e.g. a third order 2D polynomial has 10 terms)

    Returns
    -------
    poly_array : 2D numpy.ndarray

    Raises
    ------
    RuntimeError
        if the number of coefficients does not match up with a 2D polynomial
    """
    x_array = mesh_grid[0].astype('float64')
    y_array = mesh_grid[1].astype('float64')

    value = len(coeffs)
    deg = 0.0
    while value > 0.0:        
        deg += 1.0
        value -= deg
    if value < 0.0:
        raise RuntimeError('Please set number of coefficients to the factorial'
                           + 'of the polynomial degree')
    deg = int(deg - 1.0)

    poly_array = np.zeros_like(x_array)
    index = 0
    for k in xrange(deg+1):
        for j in xrange(k+1):
            i = k - j
            poly_array += coeffs[index] * x_array**i * y_array**j
            index += 1

    return poly_array

#=============================================================================#
#===== Fitting Methods =======================================================#
#=============================================================================#

def polynomial_fit(image, deg=6, center=None, radius=None, full_output=False):
    """Finds an optimal polynomial fit for an image

    Uses np.linalg.lstsq

    Args
    ----
    image : 2D numpy.ndarray
    deg : int (default=6)
        The degree of polynomial to fit.
    center : Pixel, optional
        pixel to use as the center when fitting the polynomial
    radius : number, optional
        radius from the center within which the polynomial is fit
    full_output : bool, optional
        whether or not to include the coefficients of the polynomial fit in the
        output

    Returns
    -------
    poly_fit: 2D numpy array
    coeffs : tuple
        if full_output is True
    """
    mesh_grid = mesh_grid_from_array(image)
    x_array = mesh_grid[0]
    y_array = mesh_grid[1]

    if center is not None:
        if radius is None:
            height, width = x_array.shape
            radius = min(center.x, center.y, width-center.x, height-center.y)
        x_flat = intensity_array(x_array, center, radius)
        y_flat = intensity_array(y_array, center, radius)
        image_flat = intensity_array(image, center, radius)
    else:
        x_flat = x_array.flatten()
        y_flat = y_array.flatten()
        image_flat = image.flatten()

    poly_flat = []
    for k in xrange(deg+1):
        for j in xrange(k+1):
            i = k - j
            poly_flat.append(x_flat**i * y_flat**j)
    poly_flat = np.array(poly_flat).T

    coeffs, _, _, _ = np.linalg.lstsq(poly_flat, image_flat)
    poly_fit = polynomial_array(mesh_grid, *coeffs)
    
    if center is not None:
        poly_fit = isolate_circle(poly_fit, center, radius)

    if full_output:
        return poly_fit, coeffs
    return poly_fit

def gaussian_fit(image, initial_guess=None, full_output=False, center=None, radius=None):
    """Finds an optimal gaussian fit for an image

    Uses scipy.optimize.curve_fit

    Args
    ----
    image : 2D numpy.ndarray
    initial_guess : tuple, optional
        Specifically: (x0, y0, radius, amplitude, offset)

    Returns
    -------
    gauss_fit: 2D numpy array

    """
    mesh_grid = mesh_grid_from_array(image)
    x_array = mesh_grid[0]
    y_array = mesh_grid[1]

    if center is not None and radius is not None:
        x_flat = intensity_array(x_array, center, radius)
        y_flat = intensity_array(y_array, center, radius)
        image_flat = intensity_array(image, center, radius)
    else:
        x_flat = x_array.flatten()
        y_flat = y_array.flatten()
        image_flat = image.flatten()

    if initial_guess is None:
        height, width = image.shape
        initial_guess = (width / 2.0, height / 2.0,
                         min(height, width) / 4.0,
                         image.max(),
                         image.min())

    coeffs, _ = opt.curve_fit(gaussian_array, (x_flat, y_flat),
                             image_flat, p0=initial_guess)
    gauss_fit = gaussian_array(mesh_grid, *coeffs).reshape(*image.shape)

    if center is not None and radius is not None:
        gauss_fit = isolate_circle(gauss_fit, center, radius)

    if full_output:
        return gauss_fit, coeffs
    return gauss_fit

def rectangle_fit(image, initial_guess=None, full_output=False):
    """Finds an optimal rectangle fit for an image

    Uses scipy.optimize.curve_fit

    Args
    ----
    image : 2D numpy.ndarray
    initial_guess : tuple, optional
        Specifically: (x0, y0, height, width, angle)

    Returns
    -------
    rectangle_fit: 2D numpy array
    """
    mesh_grid = mesh_grid_from_array(image)
    height, width = image.shape

    if initial_guess is None:
        initial_guess = (width / 2.0, height / 2.0,
                         width / 2.0, height / 2.0,
                         0.0)
    
    opt_parameters, _ = opt.curve_fit(rectangle_array,
                                      mesh_grid,
                                      image.ravel(),
                                      p0=initial_guess)

    rectangle_fit = rectangle_array(mesh_grid, *opt_parameters).reshape(height, width)
    if full_output:
        return rectangle_fit, opt_parameters
    return rectangle_fit
