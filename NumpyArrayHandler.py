import numpy as np
import re
from PIL import Image
from astropy.io import fits
from scipy import optimize as opt
from collections import Iterable
import matplotlib.pyplot as plt
plt.rc('font', size=16, family='sans-serif')
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)
plt.rc('lines', lw=2)

class NumpyArrayHandler(object):

    def __init__(self):
        return

    def convertImageToArray(self, image_input, full_output=False):
        """Converts an image input to a numpy array or 0.0

        Args:
            image_input: choice of (list, tuple, or individual) and (file name
                strings or numpy arrays)

        Returns:
            image_array: numpy array if types checks out or None
            bit_depth: adc conversion in bit depth
            pixel_size: length of pixel side in microns
            exp_time: length of exposure in seconds
            camera: string of the camera name
        """
        if image_input is None:
            return None

        elif isinstance(image_input, basestring):
            return self.getImageArrayFromFile(image_input, full_output)

        elif isinstance(image_input, Iterable) and isinstance(image_input[0], basestring):
            list_len = float(len(image_input))
            image_array, pixel_size, bit_depth = self.getImageArrayFromFile(image_input[0], full_output=True)
            image_array /= list_len
            for image_string in image_input[1:]:
                image_array += self.getImageArrayFromFile(image_string, full_output=False) / list_len

            if full_output:
                return image_array, pixel_size, bit_depth
            return image_array

        elif isinstance(image_input, np.ndarray) and len(image_input.shape) == 2:
            if full_output:
                return image_input, None, None
            return image_input

        elif isinstance(image_input, Iterable) and isinstance(image_input[0], np.ndarray):
            list_len = float(len(image_input))
            image_array = image_input[0] / list_len
            for image in image_input[1:]:
                image_array += image / list_len
            if full_output:
                return image_array, None, None
            return image_array
            
        else:
            raise RuntimeError('Incorrect type for image input')

    def getImageArrayFromFile(self, image_string, full_output=False):
        if image_string[-3:] == 'fit':
            image = fits.open(image_string)[0]
            image_array = image.data.astype('float64')
            header = dict(image.header)            
            bit_depth = int(header['BITPIX'])

        elif image_string[-3:] == 'tif':
            image = Image.open(image_string)
            image_array = np.array(image).astype('float64')
            bit_depth = int(image.tag[258][0])
            # Complicated way to get the header from a TIF image as a dictionary
            header = dict([i.split('=') for i in image.tag[270][0].split('\r\n')][:-1])           

        else:
            raise ValueError('Incorrect image file extension')

        if full_output:
            pixel_size = float(header['XPIXSZ'])
            #exp_time = float(header['EXPTIME'])
            #camera = str(image.header['INSTRUME'].split(':')[0])
            return image_array, pixel_size, bit_depth
        return image_array

    def getArraySum(self, image_array):
        return np.sum(image_array)

    def getRowSum(self, image_array):
        row_sum = np.sum(image_array, axis=0)
        return ((row_sum - np.min(row_sum)) / image_array.shape[0]).astype('float64')

    def getColumnSum(self, image_array):
        column_sum = np.sum(image_array, axis=1)
        return ((column_sum - np.min(column_sum)) / image_array.shape[1]).astype('float64')

    def getMeshGrid(self, image_array):
         return np.meshgrid(np.arange(image_array.shape[1]).astype('float64'),
                            np.arange(image_array.shape[0]).astype('float64'))

    def getIntensityArray(self, image_array, x0, y0, radius):
        """Finds intensities inside a circle

        Returns an array of intensities from image_array which are contained 
        within the circle centered at (x0, y0) with radius radius

        Returns:
            intensity_array: one-dimensional numpy array of intensities
        """
        image_crop, x0, y0 = self.cropImage(image_array, x0, y0, radius)
        height, width = image_crop.shape

        intensity_list = []
        for x in xrange(width):
            for y in xrange(height):
                if (x0-x)**2 + (y0-y)**2 <= (radius)**2:
                    intensity_list.append(image_crop[y,x])

        return np.array(intensity_list)

    def cropImage(self, image_array, x0, y0, radius):
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

    def removeCircle(self, image_array, x0, y0, radius, res=1):
        """Removes a circle from an array

        Args:
            radius: circle's radius
            x: horizontal pixel number for circle center
            y: vertical pixel number for circle center

        Returns:
            output_array: a copy of image_array with the circle removed
        """
        mesh_grid = self.getMeshGrid(image_array)
        return image_array * (1 - self.circleArray(mesh_grid, x0, y0, radius, res))

    def isolateCircle(self, image_array, x0, y0, radius, res=1):
        """Isolates a circle in an array

        Args:
            radius: circle's radius
            x: horizontal pixel number for circle center
            y: vertical pixel number for circle center

        Returns:
            output_array: a copy of image_array with the circle removed
        """
        mesh_grid = self.getMeshGrid(image_array)
        return image_array * self.circleArray(mesh_grid, x0, y0, radius, res)

    def applyWindow(self, image_array):        
        height, width = image_array.shape
        x_array, y_array = self.getMeshGrid(image_array)
        x0 = width/2
        y0 = height/2
        r_array = np.sqrt((x_array-x0)**2 + (y_array-y0)**2) + min(height, width) / 2
        window = self.hann_poisson_window(min(height, width), r_array)
        self.showImageArray(image_array*window)
        #self.plotOverlaidCrossSections(image_array, image_array*window, height/2, width/2)
        return image_array * window

    def hann_poisson_window(self, arr_len, arr=None):
        if arr is None:
            arr = np.arange(arr_len)
        hann = 0.5 * (1 - np.cos(2*np.pi*arr / (arr_len - 1)))
        alpha = 2.0
        poisson = np.exp(-alpha/(arr_len-1) * np.abs(arr_len - 1 - 2*arr))
        return hann * poisson

    def poisson_window(self, arr_len, arr=None):
        if arr is None:
            arr = np.arange(arr_len)
        tau = (arr_len / 2) * (8.69 / 60)
        poisson = np.exp(-np.abs(arr - (arr_len-1)/2) / tau)
        return poisson

#=============================================================================#
#===== Fitting Functions =====================================================#
#=============================================================================#

    @staticmethod
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

    @staticmethod
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
        x_array = mesh_grid[0].astype('float64')
        y_array = mesh_grid[1].astype('float64')
        x0 = float(x0)
        y0 = float(y0)
        radius = float(radius)

        if res == 1:
            circle_array = ((x_array-x0)**2 + (y_array-y0)**2 <= radius**2).astype('float64')

        else:
            circle_array = ((x_array-x0)**2 + (y_array-y0)**2 < (radius - 0.7)**2).astype('float64')

            res_array = np.arange(-0.5, 0.5, 1.0 / res) + 0.5 / res
            res_mesh_x, res_mesh_y = np.meshgrid(res_array, res_array)
            res_val = 1.0 / res**2

            for x in range(int(x0-radius), int(x0+radius) + 2):
                for y in range(int(y0-radius), int(y0+radius) + 2):
                    if circle_array[y,x] < 1.0:
                        if (x-x0)**2 + (y-y0)**2 <= (radius + 0.7)**2:
                            circle_array[y, x] += res_val * ((res_mesh_x+x-x0)**2 + (res_mesh_y+y-y0)**2 <= radius**2).astype('float64').sum()
                       
        return circle_array

    @staticmethod
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

    def getPolynomialFit(self, image_array, deg=6, x0=None, y0=None):
        """Finds an optimal polynomial fit for an image

        Returns:
            polynomial_fit: 2D numpy array
        """
        mesh_grid = self.getMeshGrid(image_array)
        if x0 is None or y0 is None:
            x0 = image_array.shape[1] / 2.0
            y0 = image_array.shape[0] / 2.0
        # Place (0, 0) at fiber center
        mesh_grid[0] -= x0
        mesh_grid[1] -= y0

        initial_guess = tuple(np.ones(deg/2 + 1)) 

        opt_coeffs, cov_matrix = opt.curve_fit(self.polynomialArray,
                                               mesh_grid,
                                               image_array.ravel(),
                                               p0=initial_guess)

        return self.polynomialArray(mesh_grid, *opt_coeffs).reshape(image_array.shape)

    def getGaussianFit(self, image_array, initial_guess=None, full_output=False):
        """Finds an optimal gaussian fit for an image

        Returns:
            polynomial_fit: 2D numpy array
        """
        mesh_grid = self.getMeshGrid(image_array)
        height, width = image_array.shape

        if initial_guess is None:
            initial_guess = (width / 2.0, height / 2.0,
                             min(height, width) / 4.0,
                             image_array.max(),
                             image_array.min())

        opt_parameters, cov_matrix = opt.curve_fit(self.gaussianArray, mesh_grid,
                                                   image_array.ravel(), p0=initial_guess)

        gaussian_fit = self.gaussianArray(mesh_grid, *opt_parameters).reshape(height, width)
        if full_output:
            return gaussian_fit, opt_parameters
        return gaussian_fit

#=============================================================================#
#===== Image Plotting ========================================================#
#=============================================================================#
    
    @staticmethod
    def plotHorizontalCrossSection(image_array, row):
        row_int = int(round(row))
        plt.plot(image_array[row_int, :])
        plt.title('Horizontal Cross Section (row = %s)'%row)
        plt.xlabel('Pixel')

    @staticmethod   
    def plotVerticalCrossSection(image_array, column):
        column_int = int(round(column))
        plt.plot(image_array[:, column_int])
        plt.title('Vertical Cross Section (column = %s)'%column)
        plt.xlabel('Pixel')

    @staticmethod
    def plotCrossSections(image_array, row, column):
        plt.figure(1)
        plt.subplot(211)
        self.plotHorizontalCrossSection(image_array, row)
        plt.subplot(212)
        self.plotVerticalCrossSection(image_array, column)
        plt.show()

    @staticmethod
    def plotOverlaidCrossSections(first_array, second_array, row, column):
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

    @staticmethod
    def plotCrossSectionSums(image_array):
        plt.figure(1)
        plt.subplot(211)
        plt.plot(self.getRowSum(image_array))
        plt.title('Average for each Column')
        plt.xlabel('Column')
        plt.subplot(212)
        plt.plot(self.getColumnSum(image_array))
        plt.title('Average for each Row')
        plt.xlabel('Row')
        plt.show()

    @staticmethod
    def showImageArray(image_array):
        plt.figure(1)
        plt.imshow(image_array, cmap='gray')
        plt.colorbar(label='intensity')
        plt.xlabel('x pixel')
        plt.ylabel('y pixel')
        plt.show()

    @staticmethod
    def show1DArray(array):
        plt.figure(1)
        plt.plot(array)
        plt.show()

    @staticmethod
    def plotFFT(freq_arrays, fft_arrays, labels=['No label']):
        plt.figure(1)
        for i in xrange(len(freq_arrays)):
            plt.plot(freq_arrays[i], fft_arrays[i], label=labels[i])
        plt.xlim(0, 1.0)
        #plt.ylim(ymax=10**-1)
        plt.yscale('log')
        plt.ylabel('Normalized Power')
        plt.xlabel('Frequency [1/um]')
        plt.legend()
        plt.show()
