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

    def gaussianArray(self, (x,y), x0, y0, radius, amp):

        gaussian_array = amp * np.exp(-2*(x - x0)**2 / radius**2
                                      -2*(y - y0)**2 / radius**2)
        gaussian_array = np.zeros([image_height, image_width])
        for x in xrange(image_width):
            for y in xrange(image_height):
                gaussian_array[y,x] = amp * np.exp(-(x - x0)**2 / (2 * (radius / 2)**2)
                                                   -(y - y0)**2 / (2 * (radius / 2)**2))
        return gaussian_array

    def 

#=============================================================================#
#===== Image Plotting ========================================================#
#=============================================================================#

    def plotHorizontalCrossSection(self, image_array, row):
        plt.plot(image_array[row, :])
        plt.title('Horizontal Cross Section (row = %s)'%row)
        plt.xlabel('Pixel')
        plt.show()

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
    