"""NumpyArrayHandler.py was written by Ryan Petersburg for use with fiber
characterization on the EXtreme PRecision Spectrograph

The functions in this module are used to handle two dimensional np.ndarray
objects that represent images from optical fibers. Functions include
converting image files to np.ndarrays, image cropping, function fitting, and 
image plotting
"""
import numpy as np
import re
from PIL import Image
from astropy.io import fits
from scipy import optimize as opt
from collections import Iterable
from datetime import datetime
from scipy.signal import medfilt2d
import matplotlib.pyplot as plt
plt.rc('font', size=20, family='serif')
plt.rc('figure', figsize=[20, 12.36])
plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('lines', lw=4)

def saveArray(input_array, file):
    """Saves a np.ndarry as the designated file

    Args:
        input_array [np.ndarray]
        file [string]
    """
    plt.imsave(file, input_array, cmap='gray')

def convertImageToArray(image_input, full_output=False):
    """Converts an image input to a numpy array or 0.0

    Args:
        image_input [list(np.ndarry or string), tuple(np.ndarray or string),
            np.ndarray, or string]: the input to be converted
        full_output [boolean]: (optional) whether or not to include relevant
            information from the image header in the return. Automatically
            False if the image_input is an np.ndarray

    Returns:
        image_array [np.ndarry or None]: 2D numpy array if the image input
            check out; None otherwise
        output_dict [dict]: dictionary of information taken from the image header
    """
    image_array = None
    output_dict = {}

    if image_input is None:
        pass

    elif isinstance(image_input, basestring):
        image_array = imageArrayFromFile(image_input, full_output)
        if full_output:
            output_dict = image_array[1]
            image_array = image_array[0]
            output_dict['num_images'] = 1

    elif isinstance(image_input, Iterable) and isinstance(image_input[0], basestring):
        list_len = float(len(image_input))
        image_array = imageArrayFromFile(image_input[0], full_output)
        if full_output:
            output_dict = image_array[1]
            image_array = image_array[0]
            output_dict['num_images'] = list_len
        image_array /= list_len
        for image_string in image_input[1:]:
            image_array += imageArrayFromFile(image_string) / list_len

    elif isinstance(image_input, Iterable) and len(np.array(image_input).shape) == 2:
        image_array = np.array(image_input)
        if full_output:
            output_dict['num_images'] = 1

    elif isinstance(image_input, Iterable) and isinstance(image_input[0], Iterable):
        list_len = float(len(image_input))
        image_input = np.array(image_input)
        image_array = image_input[0] / list_len
        for image in image_input[1:]:
            image_array += image / list_len            
        if full_output:
            output_dict['num_images'] = list_len
        
    else:
        raise RuntimeError('Incorrect type for image input')

    if full_output:
        return image_array, output_dict
    return image_array

def imageArrayFromFile(image_string, full_output=False):
    """Returns image from file as 2D np.ndarray
    
    Args:
        image_string [string]: file location to be converted
        full_output [boolean]: (optional) whether or not to include relevant
            information from the image header in the return

    Returns:
        image_array [np.ndarray]: 2D numpy array of the file's image
        output_dict [dict]: dictionary of the information from the image header
    """
    if image_string[-3:] == 'fit':
        image = fits.open(image_string)[0]
        image_array = image.data.astype('float64')
        if full_output: 
            header = dict(image.header)                     
            bit_depth = int(header['BITPIX'])

    elif image_string[-3:] == 'tif':
        image = Image.open(image_string)
        image_array = np.array(image).astype('float64')
        if full_output:
            bit_depth = int(image.tag[258][0])
            # Complicated way to get the header from a TIF image as a dictionary
            header = dict([i.split('=') for i in image.tag[270][0].split('\r\n')][:-1])           

    else:
        raise ValueError('Incorrect image file extension')

    if full_output:
        output_dict = {}
        output_dict['bit_depth'] = bit_depth
        output_dict['pixel_size'] = float(header['XPIXSZ'])
        output_dict['exp_time'] = float(header['EXPTIME'])
        #output_dict['date_time'] = datetime.strptime(header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')
        output_dict['temp'] = float(header['CCD-TEMP'])
        if 'TELESCOP' in header:
            output_dict['camera'] = str(header['TELESCOP'])
        if 'OBJECT' in header:
            output_dict['test'] = str(header['OBJECT'])

        return image_array, output_dict

    return image_array

def sumArray(image_array):
    """Sums all elements in a np.ndarray
    
    Args:
        image_array [np.ndarray]: 2D numpy array

    Returns:
        sum [float]: the sum    
    """
    return np.sum(image_array)

def sumRows(image_array):
    """Sums the rows of a 2D np.ndarray
    
    Args:
        image_array [np.ndarray]: 2D numpy array

    Returns:
        summed_row [np.ndarray]: 1D numpy array of the summed rows
    """
    row_sum = np.sum(image_array, axis=0)
    return ((row_sum - np.min(row_sum)) / image_array.shape[0]).astype('float64')

def sumColumns(image_array):
    """Sums the columnns of a 2D np.ndarray
    
    Args:
        image_array [np.ndarray]: 2D numpy array

    Returns:
        summed_column [np.ndarray]: 1D numpy array of the summed columns
    """
    column_sum = np.sum(image_array, axis=1)
    return ((column_sum - np.min(column_sum)) / image_array.shape[1]).astype('float64')

def meshGridFromArray(image_array):
    """Creates a numpy meshgrid of pixel number for an image

    Args:
        image_array: 2D numpy image

    Returns:
        mesh_grid: (2D array of x values, 2D array of y values)
    """
    return np.meshgrid(np.arange(image_array.shape[1]).astype('float64'),
                       np.arange(image_array.shape[0]).astype('float64'))

def intensityArray(image_array, x0, y0, radius):
    """Finds intensities inside a circle

    Returns an array of intensities from image_array which are contained 
    within the circle centered at (x0, y0) with radius radius

    Returns:
        intensity_array: one-dimensional numpy array of intensities
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

    Returns:
        image_crop: cropped image
        y0: new center y
        x0: new center x
    """
    image_crop = image_array[int(y0-radius):int(y0+radius)+2,
                             int(x0-radius):int(x0+radius)+2]
    y0 = radius + (y0-radius)-int(y0-radius)
    x0 = radius + (x0-radius)-int(x0-radius)
    return image_crop, x0, y0

def removeCircle(image_array, x0, y0, radius, res=1):
    """Removes a circle from an array

    Args:
        radius: circle's radius in pixels
        x: horizontal pixel number for circle center
        y: vertical pixel number for circle center

    Returns:
        output_array: a copy of image_array with the circle removed
    """
    mesh_grid = meshGridFromArray(image_array)
    return image_array * (1 - circleArray(mesh_grid, x0, y0, radius, res))

def isolateCircle(image_array, x0, y0, radius, res=1):
    """Isolates a circle in an array

    Args:
        radius: circle's radius in pixels
        x: horizontal pixel number for circle center
        y: vertical pixel number for circle center

    Returns:
        output_array: a copy of image_array with the circle isolated
    """
    mesh_grid = meshGridFromArray(image_array)
    return image_array * circleArray(mesh_grid, x0, y0, radius, res)

def applyWindow(image_array):
    """Applies a FFT window to an image
    
    Args:
        image_array [np.ndarray]: 2D numpy array

    Returns:
        windowed_array [np.ndarray]: 2D numpy array of the windowed image_array
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

    Args:
        arr_len [int]: length of the array on which the window is built
        arr [np.ndarray]: (optional) any dimensional numpy array on which
            the window is applied

    Returns:
        hann_poisson_window [np.ndarray]: 1D numpy array
    """
    if arr is None:
        arr = np.arange(arr_len)
    hann = 0.5 * (1 - np.cos(2*np.pi*arr / (arr_len - 1)))
    alpha = 2.0
    poisson = np.exp(-alpha/(arr_len-1) * np.abs(arr_len - 1 - 2*arr))
    return hann * poisson

def poisson_window(arr_len, arr=None):
    if arr is None:
        arr = np.arange(arr_len)
    tau = (arr_len / 2) * (8.69 / 60)
    poisson = np.exp(-np.abs(arr - (arr_len-1)/2) / tau)
    return poisson

def filteredImage(image_array, kernel_size):
    """Applies a median filter to an image
    
    Args:
        image_array [np.ndarray]: 2D numpy array
        kernel_size [int]: odd number that gives the side length of the kernel
            which the median filter uses

    Returns:
        filtered_image [np.ndarray]: 2D numpy array of the filtered image
    """
    if kernel_size < 2.0:
        return image_array
    return medfilt2d(image_array, kernel_size)

#=============================================================================#
#===== 2D Array Functions ====================================================#
#=============================================================================#

def gaussianArray(mesh_grid, x0, y0, radius, amp, offset):
    """Creates a 2D gaussian function as a 1D array

    Args:
        mesh_grid: independent variables x and y separated into two arrays
            each with the proper dimensions (np.meshgrid)
        x0: center position x of gaussian
        y0: center position y of gaussian
        radius: radius of gaussian (2 standard deviations or 95% volume)
        amp: amplitude of gaussian

    Returns:
        Ravelled gaussian numpy array (single dimension) usable in
            Scipy.optimize.curve_fit method and properly reshaped by
            gaussian_array.reshape(height, width)
    """
    gaussian_array = offset + amp * np.exp(-2*(mesh_grid[0] - x0)**2 / radius**2
                                           -2*(mesh_grid[1] - y0)**2 / radius**2)
    return gaussian_array.ravel()

def circleArray(mesh_grid, x0, y0, radius, res=1):
    """Creates a 2D tophat function of amplitude 1.0
    
    Args:
        res: Resolution element for more precise circular edges


    Returns:
        circle_array: 2D numpy array of float values where points inside
            the circle are 1.0 and outside the circle are 0.0. Points along
            the edge of the circle are weighted based on their relative
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

        for x in range(int(x0-radius), int(x0+radius) + 2):
            for y in range(int(y0-radius), int(y0+radius) + 2):
                if circle_array[y, x] < 1.0:
                    if (x-x0)**2 + (y-y0)**2 <= (radius + np.sqrt(2) / 2.0)**2:
                        circle_array[y, x] += res_val * ((res_mesh_x+x-x0)**2 + (res_mesh_y+y-y0)**2 <= radius**2).astype('float64').sum()
    return circle_array

def polynomialArray(mesh_grid, *coeff):
    """Creates an even 2D polynomial of arbitrary degree for given x, y

    Uses a mesh grid and list of coefficients to create a two dimensional
    even polynomial array that is azimuthally symmetric around the (0, 0)
    point of the mesh grid, e.g. 3*r^4 + 2*r^2 + 1 where
    r = sqrt(x^2 + y^2)

    Args:
        *coeffs: coefficients where the length is the degree divided by two
            (since the polynomial is even)

    Returns:
        polynomial_array: 2D numpy array of float values
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

    Returns:
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

    opt_parameters, cov_matrix = opt.curve_fit(gaussianArray, mesh_grid,
                                               image_array.ravel(), p0=initial_guess)

    gaussian_fit = gaussianArray(mesh_grid, *opt_parameters).reshape(height, width)
    if full_output:
        return gaussian_fit, opt_parameters
    return gaussian_fit

#=============================================================================#
#===== Image Plotting ========================================================#
#=============================================================================#

def plotHorizontalCrossSection(image_array, row):
    row_int = int(round(row))
    plt.plot(image_array[row_int, :])
    plt.title('Horizontal Cross Section (row = %s)'%row)
    plt.xlabel('Pixel')
  
def plotVerticalCrossSection(image_array, column):
    column_int = int(round(column))
    plt.plot(image_array[:, column_int])
    plt.title('Vertical Cross Section (column = %s)'%column)
    plt.xlabel('Pixel')

def plotCrossSections(image_array, row, column):
    plt.figure(1)
    plt.subplot(211)
    NumpyArrayHandler.plotHorizontalCrossSection(image_array, row)
    plt.subplot(212)
    NumpyArrayHandler.plotVerticalCrossSection(image_array, column)
    plt.show()

def plotOverlaidCrossSections(first_array, second_array, row, column):
    row = int(round(row))
    column = int(round(column))
    plt.figure(1)
    plt.subplot(211)
    plt.plot(first_array[row, :])
    plt.plot(second_array[row, :])
    plt.title('Horizontal Cross Section (row = %s)'%row)
    plt.xlabel('Pixel')
    plt.subplot(212)
    plt.plot(first_array[:, column])
    plt.plot(second_array[:, column])
    plt.title('Vertical Cross Section (column = %s)'%column,)
    plt.xlabel('Pixel')
    plt.show()

def plotCrossSectionSums(image_array):
    plt.figure(1)
    plt.subplot(211)
    plt.plot(sumRows(image_array))
    plt.title('Average for each Column')
    plt.xlabel('Column')
    plt.subplot(212)
    plt.plot(sumColumns(image_array))
    plt.title('Average for each Row')
    plt.xlabel('Row')
    plt.show()

def showImageArray(image_array):
    plt.figure(1)
    plt.imshow(image_array, cmap='gray')
    plt.colorbar(label='intensity')
    plt.xlabel('x pixel')
    plt.ylabel('y pixel')
    plt.show()

def show1DArray(array):
    plt.figure(1)
    plt.plot(array)
    plt.show()

def plotFFT(freq_arrays, fft_arrays, labels=['No label'], title='Power Spectrum'):
    plt.figure(1)
    for i in xrange(len(freq_arrays)):
        plt.plot(freq_arrays[i], fft_arrays[i], label=labels[i])
    plt.xlim(0, freq_arrays[0].max()/2.0)
    plt.yscale('log')
    plt.ylabel('Normalized Power')
    plt.xlabel('Frequency [1/um]')
    plt.title(title)
    plt.legend()
    plt.show()
