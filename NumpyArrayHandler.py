import numpy as np
from PIL import Image
from scipy import optimize as opt
import itertools
import matplotlib.pyplot as plt
import time

class NumpyArrayHandler(object):

    def __init__(self):
        return

    def convertImageToArray(self, image_string):
        return np.array(Image.open(image_string))[:, :, 0].astype(float)    

    def getArraySum(self, image_array):
        return np.sum(image_array)

    def getRowSum(self, image_array):
        row_sum = np.sum(image_array, axis=0)
        return ((row_sum - np.min(row_sum)) / image_array.shape[0]).astype(float)

    def getColumnSum(self, image_array):
        column_sum = np.sum(image_array, axis=1)
        return ((column_sum - np.min(column_sum)) / image_array.shape[1]).astype(float)

    def getMeshGrid(self, image_array):
         return np.meshgrid(np.arange(image_array.shape[1]).astype('float64'),
                            np.arange(image_array.shape[0]).astype('float64'))

#=============================================================================#
#===== Fitting Functions =====================================================#
#=============================================================================#

    def gaussianArray(self, mesh_grid, x0, y0, radius, amp):
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
        gaussian_array = amp * np.exp(-2*(mesh_grid[1] - x0)**2 / radius**2
                                      -2*(mesh_grid[0] - y0)**2 / radius**2)
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
            circle_array = 1.0 * ((x_array-x0)**2 + (y_array-y0)**2 <= radius**2).astype(float)

        else:
            circle_array = 1.0 * ((x_array-x0)**2 + (y_array-y0)**2 < (radius - 0.7)**2).astype(float)

            res_array = np.arange(-0.5, 0.5, 1.0 / res) + 0.5 / res
            res_mesh_x, res_mesh_y = np.meshgrid(res_array, res_array)
            res_val = 1.0 / res**2

            for x in range(int(x0-radius), int(x0+radius) + 2):
                for y in range(int(y0-radius), int(y0+radius) + 2):
                    if circle_array[y,x] == 0.0:
                        if (x-x0)**2 + (y-y0)**2 <= (radius + 0.7)**2:
                            circle_array[y, x] = res_val * ((res_mesh_x+x-x0)**2 + (res_mesh_y+y-y0)**2 <= radius**2).astype(float).sum()
                       
        return circle_array

    def polynomialArray(self, mesh_grid, *coeff):
        """Creates a 2D polynomial function of arbitrary degree

        Args:
            *coeffs: coefficients where the number determines the degree

        Returns:
            polynomial_array: 2D numpy array of float values
        """
        x_array = mesh_grid[0]
        y_array = mesh_grid[1]

        polynomial_array = np.zeros_like(x_array)
        deg = int(np.sqrt(len(coeff)) - 1)

        print coeff

        ij = itertools.product(range(deg+1), range(deg+1))
        for d, (i,j) in enumerate(ij):
            polynomial_array += coeff[d] * (x_array)**i * (y_array)**j

        return polynomial_array.ravel()

#=============================================================================#
#===== Fitting Methods =======================================================#
#=============================================================================#

    def getPolynomialFit(self, image_array, deg=4, x0=None, y0=None):
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

        initial_guess = (image_array.max(),) + tuple(-1.0 * np.ones((deg+1)**2 - 1))        

        opt_coeffs, cov_matrix = opt.curve_fit(self.polynomialArray,
                                               mesh_grid,
                                               image_array.ravel(),
                                               p0=initial_guess)

        return self.polynomialArray(mesh_grid, *opt_coeffs).reshape(image_array.shape)

    def getPolynomialFit2(self, image_array, deg=4):
        height, width = image_array.shape
        mesh_grid = self.getMeshGrid(image_array)
        x_array = mesh_grid[0].astype('float64')
        y_array = mesh_grid[1].astype('float64')
        num_terms = (deg+1)**2

        poly_tensor = np.zeros((height, width, num_terms))
        ij = itertools.product(range(deg+1), range(deg+1))
        for k, (i,j) in enumerate(ij):
            poly_tensor[:,:,k] = x_array**i * y_array**j

        opt_coeffs, resids, rk, s = np.linalg.lstsq(poly_tensor.reshape(num_terms, height*width).swapaxes(0,1),
                                                    image_array.reshape(height*width))

        return self.polynomialArray(mesh_grid, *opt_coeffs).reshape(image_array.shape)

    def getGaussianFit(self, image_array, initial_guess=None, full_output=False):
        mesh_grid = self.getMeshGrid(image_array)
        height, width = image_array.shape

        if initial_guess is None:
            initial_guess = (width / 2.0, height / 2.0,
                             min(height, width) / 4.0, image_array.max())

        opt_parameters, cov_matrix = opt.curve_fit(self.gaussianArray, mesh_grid,
                                                   image_array.ravel(), p0=initial_guess)

        gaussian_fit = self.gaussianArray(self.mesh_grid, *opt_parameters).reshape(self.image_height,
                                                                                   self.image_width)
        if full_output:
            return gaussian_fit, opt_parameters
        return gaussian_fit

#=============================================================================#
#===== Image Plotting ========================================================#
#=============================================================================#

    def plotHorizontalCrossSection(self, image_array, row):
        plt.plot(image_array[row, :])
        plt.title('Horizontal Cross Section (row = %s)'%row)
        plt.xlabel('Pixel')

    def plotVerticalCrossSection(self, image_array, column):
        plt.plot(image_array[:, column])
        plt.title('Vertical Cross Section (column = %s)'%column)
        plt.xlabel('Pixel')

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
        plt.plot(first_array[row, :])
        plt.plot(second_array[row, :])
        plt.title('Horizontal Cross Section (row = %s)'%row)
        plt.xlabel('Pixel')
        plt.subplot(212)
        plt.plot(first_array[:, column])
        plt.plot(second_array[:, column])
        plt.title('Vertical Cross Section (column = %s)'%column)
        plt.xlabel('Pixel')
        plt.show()

    def plotCrossSectionSums(self, image_array):
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

    def showImageArray(self, image_array):
        plt.figure(1)
        plt.imshow(image_array, cmap='gray')
        plt.colorbar(label='intensity')
        plt.xlabel('x pixel')
        plt.ylabel('y pixel')
        plt.show()

    def show1DArray(self, array):
        plt.figure(1)
        plt.plot(array)
        plt.show()
