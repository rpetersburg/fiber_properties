"""NumpyArrayHandler.py was written by Ryan Petersburg for use with fiber
characterization on the EXtreme PREcision Spectrograph

The functions in this module are used to handle two dimensional np.ndarray
objects that represent images from optical fibers. Functions include
array summing, image cropping, image filtering, function fitting, and
function generation.
"""
import numpy as np
from scipy import optimize as opt
from scipy.signal import medfilt2d
from astropy.io import fits

#=============================================================================#
#===== Array Summing =========================================================#
#=============================================================================#

def sumArray(image_array):
    """Returns the sum of all elements in a numpy.ndarray"""
    return np.sum(image_array)

def sumRows(image_array):
    """Sums the rows of a 2D np.ndarray
    
    Args
    ----
    image_array : 2D numpy.ndarray

    Returns
    -------
    summed_row : 1D numpy.ndarray

    """
    row_sum = np.sum(image_array, axis=0)
    return ((row_sum - np.min(row_sum)) / image_array.shape[0]).astype('float64')

def sumColumns(image_array):
    """Sums the columnns of a 2D np.ndarray
    
    Args
    ----
    image_array : 2D numpy.ndarray

    Returns
    -------
    summed_column : 1D numpy.ndarray
    """
    column_sum = np.sum(image_array, axis=1)
    return ((column_sum - np.min(column_sum)) / image_array.shape[1]).astype('float64')

#=============================================================================#
#===== Array Alterations =====================================================#
#=============================================================================#

def meshGridFromArray(image_array):
    """Creates a numpy meshgrid of pixel number for an image

    Args
    ----
    image_array : 2D numpy.ndarray

    Returns
    -------
    mesh_grid_x : 2D numpy.ndarray
        The x position of each point in the grid
    mesh_grid_y : 2D numpy.ndarray
        The y position of each point in the grid
    """
    return np.meshgrid(np.arange(image_array.shape[1]).astype('float64'),
                       np.arange(image_array.shape[0]).astype('float64'))

def intensityArray(image_array, x0, y0, radius):
    """Returns intensities from inside a circle

    Returns an array of intensities from image_array which are contained 
    within the circle with radius centered at (x0, y0)

    Args
    ----
    image_array : 2D numpy.ndarray
    x0 : number (pixels)
    y0 : number (pixels)
    radius : number (pixels)

    Returns
    -------
    intensity_array : 1D numpy.ndarray
        Intensities of the elements contained within the given circle
    """
    image_crop, x0, y0 = cropImage(image_array, x0, y0, radius)
    height, width = image_crop.shape

    intensity_list = []
    for x in xrange(width):
        for y in xrange(height):
            if (x0-x)**2 + (y0-y)**2 <= (radius)**2:
                intensity_list.append(image_crop[y,x])

    return np.array(intensity_list)

def cropImage(image_array, x0, y0, radius):
    """Crops image to square with radius centered at (y0, x0)

    Args
    ----
    image_array : 2D numpy.ndarray
    x0 : number (pixels)
    y0 : number (pixels)
    radius : number (pixels)

    Returns
    -------
    image_crop : 2D numpy.ndarray
    new_x0 : float
    new_y0 : float
    """
    image_crop = image_array[int(y0-radius):int(y0+radius)+2,
                             int(x0-radius):int(x0+radius)+2]
    new_y0 = y0 - int(y0-radius)
    new_x0 = x0 - int(x0-radius)

    return image_crop, new_x0, new_y0

def subframeImage(image, subframe_x, subframe_y, width, height):
    """Creates the subframe of an image with the given parameters."""
    return image[subframe_y : subframe_y + height,
                 subframe_x : subframe_x + width]

def removeCircle(image_array, x0, y0, radius, res=1):
    """Removes a circle from an array

    Args
    ----
    image_array : 2D numpy.ndarray
    x0 : number (pixels)
    y0 : number (pixels)
    radius : number (pixels)
    
    Returns
    -------
    removed_circle_array : 2D numpy.ndarray
        Input image array with the defined circle removed
    """
    mesh_grid = meshGridFromArray(image_array)
    return image_array * (1 - circleArray(mesh_grid, x0, y0, radius, res))

def isolateCircle(image_array, x0, y0, radius, res=1):
    """Isolates a circle in an array

    Args
    ----
    image_array : 2D numpy.ndarray
    x0 : number (pixels)
    y0 : number (pixels)
    radius : number (pixels)

    Returns
    -------
    isolated_circle_array : 2D numpy.ndarray
        Input image array with the defined circle isolated in the image
    """
    mesh_grid = meshGridFromArray(image_array)
    return image_array * circleArray(mesh_grid, x0, y0, radius, res)

def applyWindow(image_array):
    """Applies a FFT window to an image
    
    Args
    ----
    image_array : 2D numpy.ndarray

    Returns
    -------
    windowed_array : 2D numpy.ndarray
    """        
    height, width = image_array.shape
    x_array, y_array = meshGridFromArray(image_array)
    x0 = width/2
    y0 = height/2
    r_array = np.sqrt((x_array-x0)**2 + (y_array-y0)**2) + min(height, width) / 2
    window = hann_poisson_window(min(height, width), r_array)
    return image_array * window

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

def filteredImage(image_array, kernel_size):
    """Applies a median filter to an image
    
    Args
    ----
    image_array : 2D numpy.ndarray
    kernel_size : int (odd)
        side length of the kernel which the median filter uses

    Returns
    -------
    filtered_image : 2D numpy.ndarray

    """
    if kernel_size < 2.0:
        return image_array
    return medfilt2d(image_array, kernel_size)

#=============================================================================#
#===== 2D Array Functions ====================================================#
#=============================================================================#

def gaussianArray(mesh_grid, x0, y0, radius, amp, offset):
    """Creates a 2D gaussian function as a 1D array

    Args
    ----
    mesh_grid : numpy.meshgrid
    image_array : 2D numpy.ndarray
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

def circleArray(mesh_grid, x0, y0, radius, res=1):
    """Creates a 2D tophat function of amplitude 1.0
    
    Args
    ----
    mesh_grid : numpy.meshgrid
    image_array : 2D numpy.ndarray
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
        circle_array = ((x_array-x0)**2 + (y_array-y0)**2 < (radius - np.sqrt(2) / 2.0)**2).astype('float64')

        res_array = np.arange(-0.5, 0.5, 1.0 / res) + 0.5 / res
        res_mesh_x, res_mesh_y = np.meshgrid(res_array, res_array)
        res_val = 1.0 / res**2
        height, width = circle_array.shape

        for x in range(max(0, int(x0-radius)), min(width, int(x0+radius) + 2)):
            for y in range(max(0, int(y0-radius)), min(height, int(y0+radius) + 2)):
                if circle_array[y, x] < 1.0:
                    if (x-x0)**2 + (y-y0)**2 <= (radius + np.sqrt(2) / 2.0)**2:
                        circle_array[y, x] += res_val * ((res_mesh_x+x-x0)**2 + (res_mesh_y+y-y0)**2 <= radius**2).astype('float64').sum()
    return circle_array

def polynomialArray(mesh_grid, *coeff):
    """Even, 2D, radial polynomial of arbitrary degree for given x, y

    Uses a mesh grid and list of coefficients to create a two dimensional
    even polynomial array that is azimuthally symmetric around the (0, 0)
    point of the mesh grid, e.g. 3*r^4 + 2*r^2 + 1 where
    r = sqrt(x^2 + y^2)

    Args
    ----
    *coeffs :
        coefficients where the length is the degree divided by two
        (since the polynomial is even)

    Returns
    -------
    polynomial_array : 2D numpy.ndarray

    """
    x_array = mesh_grid[0]
    y_array = mesh_grid[1]
    r_array = x_array**2 + y_array**2

    polynomial_array = np.zeros_like(r_array)
    for i in xrange(len(coeff)):
        polynomial_array += coeff[i] * r_array**i

    return polynomial_array.ravel()

#=============================================================================#
#===== Fitting Methods =======================================================#
#=============================================================================#

def polynomialFit(image_array, deg=6, x0=None, y0=None):
    """Finds an optimal polynomial fit for an image

    Uses scipy.optimize.curve_fit

    Args
    ----
    image_array : 2D numpy.ndarray
    deg : int (default=6)
        The degree of polynomial to fit.
    x0 : number
        The center column to use for the radial polynomial. Uses center of
        image_array if None.
    y0 : number
        The center row to use for the radial polynomial. Uses center of
        image_array if None.

    Returns
    -------
    polynomial_fit: 2D numpy array

    """
    mesh_grid = meshGridFromArray(image_array)
    if x0 is None or y0 is None:
        x0 = image_array.shape[1] / 2.0
        y0 = image_array.shape[0] / 2.0
    # Place (0, 0) at fiber center
    mesh_grid[0] -= x0
    mesh_grid[1] -= y0

    initial_guess = tuple(np.ones(deg/2 + 1)) 

    opt_coeffs, cov_matrix = opt.curve_fit(polynomialArray,
                                           mesh_grid,
                                           image_array.ravel(),
                                           p0=initial_guess)

    return polynomialArray(mesh_grid, *opt_coeffs).reshape(image_array.shape)

def gaussianFit(image_array, initial_guess=None, full_output=False):
    """Finds an optimal gaussian fit for an image

    Uses scipy.optimize.curve_fit

    Args
    ----
    image_array : 2D numpy.ndarray
    initial_guess : tuple, optional
        Specifically: (x0, y0, radius, amplitude, offset)

    Returns:
        polynomial_fit: 2D numpy array

    """
    mesh_grid = meshGridFromArray(image_array)
    height, width = image_array.shape

    if initial_guess is None:
        initial_guess = (width / 2.0, height / 2.0,
                         min(height, width) / 4.0,
                         image_array.max(),
                         image_array.min())

    opt_parameters, _ = opt.curve_fit(gaussianArray,
                                      mesh_grid,
                                      image_array.ravel(),
                                      p0=initial_guess)

    gaussian_fit = gaussianArray(mesh_grid, *opt_parameters).reshape(height, width)
    if full_output:
        return gaussian_fit, opt_parameters
    return gaussian_fit