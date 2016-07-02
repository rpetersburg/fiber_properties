import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

class NumpyArrayHandler(object):

    def __init__(self):
        return

    def convertImageToArray(self, image_string):
        return np.array(Image.open(image_string))[:, :, 0].astype(float)    

    def getArraySum(self, image_array):
        return np.sum(image_array)

    def getRowSum(self, image_array):
        row_sum = np.sum(image_array, axis=0)
        return ((row_sum - np.min(row_sum)) / self.image_height).astype(float)

    def getColumnSum(self, image_array):
        column_sum = np.sum(image_array, axis=1)
        return ((column_sum - np.min(column_sum)) / self.image_width).astype(float)

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
            Ravelled tophat numpy array (single dimension) usable in
                Scipy.optimize.curve_fit method and properly reshaped by
                gaussian_array.reshape(height, width)
        """
        gaussian_array = amp * np.exp(-2*(mesh_grid[0] - x0)**2 / radius**2
                                      -2*(mesh_grid[1] - y0)**2 / radius**2)
        return gaussian_array.ravel()

    def tophatArray(self, mesh_grid, x0, y0, radius, amp):
        """Creates a 2D tophat function as a 1D array

        Args:
            mesh_grid: independent variables x and y separated into two arrays
                each with the proper dimensions (np.meshgrid)
            x0: center position x of tophat
            y0: center position y of tophat
            radius: radius of tophat
            amp: amplitude of tophat

        Returns:
            tophat_array: Ravelled tophat numpy array (single dimension) usable in
                Scipy.optimize.curve_fit method and properly reshaped by
                tophat_array.reshape(height, width)
        """
        res = 1 # Resolution element for more precise circular edges

        x_array = mesh_grid[0]
        y_array = mesh_grid[1]
        x0 = float(x0)
        y0 = float(y0)
        radius = float(radius)
        amp = float(amp)

        if res == 1:
            tophat_array = amp * ((x_array-x0)**2 + (y_array-y0)**2 <= (radius * 1.001)**2).astype(float)
        else:
            tophat_array = amp * ((x_array-x0)**2 + (y_array-y0)**2 < (radius-1)**2).astype(float)

            height, width = tophat_array.shape
            for x in xrange(width):
                for y in xrange(height):
                    if (x-x0)**2 + (y-y0)**2 >= (radius-1)**2 and (x-x0)**2 + (y-y0)**2 <= (radius+1)**2:
                        for i in xrange(res):
                            for j in xrange(res):
                                x_temp = x - 0.5 + 1 / (2.0*res) + i / float(res)
                                y_temp = y - 0.5 + 1 / (2.0*res) + j / float(res)
                                if (x_temp-x0)**2 + (y_temp-y0)**2 < radius**2:
                                    tophat_array[y, x] += amp / res**2

        return tophat_array

    def circleArray(self, mesh_grid, x0, y0, radius):
        """Creates a 2D tophat function of amplitude 1.0

         Returns:
            circle_array: 2D numpy array of float values where points inside the circle are
                1.0 and outside the circle are 0.0. Points along the edge of the circle
                are weighted based on their relative distance to the center
        """
        res = 1 # Resolution element for more precise circular edges

        x_array = mesh_grid[0]
        y_array = mesh_grid[1]
        x0 = float(x0)
        y0 = float(y0)
        radius = float(radius)

        if res == 1:
            circle_array = 1.0 * ((x_array-x0)**2 + (y_array-y0)**2 <= (radius * 1.001)**2).astype(float)

        else:
            circle_array = 1.0 * ((x_array-x0)**2 + (y_array-y0)**2 < (radius-1)**2).astype(float)

            height, width = circle_array.shape
            for x in xrange(width):
                for y in xrange(height):
                    if (x-x0)**2 + (y-y0)**2 >= (radius-1)**2 and (x-x0)**2 + (y-y0)**2 <= (radius+1)**2:
                        for i in xrange(res):
                            for j in xrange(res):
                                x_temp = x - 0.5 + 1 / (2.0*res) + i / float(res)
                                y_temp = y - 0.5 + 1 / (2.0*res) + j / float(res)
                                if (x_temp-x0)**2 + (y_temp-y0)**2 < (radius)**2:
                                    circle_array[y, x] += 1.0 / res**2

        return circle_array      

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
        plt.imshow(image_array)
        plt.show()
    