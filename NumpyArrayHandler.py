import numpy as np

class NumpyArrayHandler(object):

    def __init__(self):

    def getArraySum(self, image_array):
        return np.sum(image_array)

    def getRowSum(self, image_array):
        row_sum = np.sum(image_array, axis=0)
        return ((row_sum - np.min(row_sum)) / self.image_height).astype(float)

    def getColumnSum(self, image_array):
        column_sum = np.sum(image_array, axis=1)
        return ((column_sum - np.min(column_sum)) / self.image_width).astype(float)

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
    