import numpy as np
from PIL import Image
from scipy import optimize as opt
import itertools
import time
import matplotlib.pyplot as plt
font = {'size' : 14,
        'family' : 'serif'}
lines = {'lw' : 2}
plt.rc('font', **font)
plt.rc('lines', **lines)

class NumpyArrayHandler(object):

    def __init__(self):
        return

    def convertImageToArray(self, image_input):
        """Converts an image input to a numpy array or 0.0

        Args:
            image_input: choice of (list, tuple, or individual) and (file name
                strings or numpy arrays)

        Returns:
            image_array: numpy array if types checks out
            None: if there's no input, return None
        """
        if not image_input:
            return None

        elif isinstance(image_input, basestring):
            return np.array(Image.open(image_input)).astype('float64')

        elif isinstance(image_input, (list, tuple)) and isinstance(image_input[0], basestring):
            image_array = np.zeros_like(np.array(Image.open(image_input[0])).astype('float64'))
            for image_string in image_input:
                image_array += np.array(Image.open(image_string)).astype('float64') / float(len(image_input))
            return image_array

        elif isinstance(image_input, np.ndarray):
            return image_input

        elif isinstance(image_input, (list, tuple)) and isinstance(image_input[0], np.ndarray):
            image_array = np.zeros_like(image_input[0])
            for image in image_input:
                image_array += image / float(len(image_input))
            return image_array
            
        else:
            raise ValueError('Incorrect type for image input')

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

#=============================================================================#
#===== Fitting Functions =====================================================#
#=============================================================================#

    def gaussianArray(self, mesh_grid, x0, y0, radius, amp, offset):
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

    def circleArray(self, mesh_grid, x0, y0, radius, res=1):
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
            circle_array = 1.0 * ((x_array-x0)**2 + (y_array-y0)**2 <= radius**2).astype('float64')

        else:
            circle_array = 1.0 * ((x_array-x0)**2 + (y_array-y0)**2 < (radius - 0.7)**2).astype('float64')

            res_array = np.arange(-0.5, 0.5, 1.0 / res) + 0.5 / res
            res_mesh_x, res_mesh_y = np.meshgrid(res_array, res_array)
            res_val = 1.0 / res**2

            for x in range(int(x0-radius), int(x0+radius) + 2):
                for y in range(int(y0-radius), int(y0+radius) + 2):
                    if circle_array[y,x] == 0.0:
                        if (x-x0)**2 + (y-y0)**2 <= (radius + 0.7)**2:
                            circle_array[y, x] = res_val * ((res_mesh_x+x-x0)**2 + (res_mesh_y+y-y0)**2 <= radius**2).astype('float64').sum()
                       
        return circle_array

    def polynomialArray(self, mesh_grid, *coeff):
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

    def plotHorizontalCrossSection(self, image_array, row):
        row_int = int(round(row))
        plt.plot(image_array[row_int, :], linewidth=4)
        plt.title('Horizontal Cross Section (row = %s)'%row, fontsize='large')
        plt.tick_params(labelsize='large')
        plt.xlabel('Pixel', fontsize='large')

    def plotVerticalCrossSection(self, image_array, column):
        column_int = int(round(column))
        plt.plot(image_array[:, column_int], linewidth=4)
        plt.title('Vertical Cross Section (column = %s)'%column, fontsize='large')
        plt.tick_params(labelsize='large')
        plt.xlabel('Pixel', fontsize='large')

    def plotCrossSections(self, image_array, row, column):
        plt.figure(1)
        plt.subplot(211)
        self.plotHorizontalCrossSection(image_array, row)
        plt.subplot(212)
        self.plotVerticalCrossSection(image_array, column)
        plt.show()

    def plotOverlaidCrossSections(self, first_array, second_array, row, column):
        plt.figure(1)
        plt.subplot(211)
        plt.plot(first_array[row, :], linewidth=4)
        plt.plot(second_array[row, :], linewidth=4)
        plt.title('Horizontal Cross Section (row = %s)'%row, fontsize='large')
        plt.xlabel('Pixel', fontsize='large')
        plt.subplot(212)
        plt.plot(first_array[:, column], linewidth=4)
        plt.plot(second_array[:, column], linewidth=4)
        plt.title('Vertical Cross Section (column = %s)'%column, fontsize='large')
        plt.xlabel('Pixel', fontsize='large')
        plt.show()

    def plotCrossSectionSums(self, image_array):
        plt.figure(1)
        plt.subplot(211)
        plt.plot(self.getRowSum(image_array), linewidth=4)
        plt.title('Average for each Column', fontsize='large')
        plt.tick_params(labelsize='large')
        plt.xlabel('Column', fontsize='large')
        plt.subplot(212)
        plt.plot(self.getColumnSum(image_array), linewidth=4)
        plt.title('Average for each Row', fontsize='large')
        plt.xlabel('Row', fontsize='large')
        plt.show()

    def showImageArray(self, image_array):
        plt.figure(1)
        plt.imshow(image_array, cmap='gray')
        plt.colorbar(label='intensity')
        plt.xlabel('x pixel', fontsize='large')
        plt.ylabel('y pixel', fontsize='large')
        plt.show()

    def show1DArray(self, array):
        plt.figure(1)
        plt.plot(array, linewidth=4)
        plt.tick_params(labelsize='large')
        plt.show()
