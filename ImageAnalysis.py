import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import NumpyArrayHandler

class ImageAnalysis(NumpyArrayHandler):

    def __init__(self, image_string, dark_image_strings=[None],
                                     flat_field_image_strings=[None],
                                     backlit_image_string=[None]):
        self.threshold = 1

        self.setImageArray(image_string)
        self.setDarkImage(*dark_image_strings)
        self.setFlatFieldImage(*flat_field_image_strings)
        self.setBacklitImage(backlit_image_string)
        """
        if dark_image_strings is not None:
            self.setDarkImage(*dark_image_strings)
        else:
            self.dark_image_array = None
        if flat_field_image_strings is not None:
            self.setFlatFieldImage(*flat_field_image_strings)
        else:
            self.flat_field_image_array = None
        if backlit_image_string is not None:
            self.setBacklitImage(backlit_image_string)
        else:
            self.backlit_image_array = None
        """

        self.image_height, self.image_width = self.image_array.shape

        self.left = None
        self.right = None
        self.top = None
        self.bottom = None
        self.fiber_diameter = None

        self.centroid_y = None
        self.centroid_x = None
        self.center_y_circle = None
        self.center_x_circle = None
        self.center_y_edge = None
        self.center_x_edge = None

        # Golden Ratio for maximization tests
        self.phi = (5 ** 0.5 - 1) / 2

#==============================================================================#
#==== Image Initialization ====================================================#
#==============================================================================#

    def convertImageToArray(self, image_string):
        return np.array(Image.open(image_string))[:, :, 0].astype(float)

    def executeErrorCorrections(self, image_array=None):
        """Applies corrective images to fiber image

        Uses Flat Field Image and Dark Image to correct for errors in the
        detector. If called with an argument, executes the corrections to the
        given image array. Otherwise, executes the corrections using the fiber
        image initialized with the ImageAnalysis instance.

        Args:
            image_array: numpy array of the fiber image or None to use the fiber
                image contained in the class instance

        Returns:
            If an argument is given, the corrected image numpy array
            If no argument is given, None
        """
        if self.flat_field_image_array is None or self.dark_image_array is None:
            if self.flat_field_image_array is None:
                print 'Flat Field Image has not been set'
            if self.dark_image_array is None:
                print 'Dark Image has not been set'
            return

        if image_array is None:
            self.image_array = self.removeDarkImage(self.image_array)
            self.image_array = self.image_array \
                               * np.mean(self.flat_field_image_array) \
                               / (self.flat_field_image_array
                                  - self.dark_image_array)
            return
        else:
            image_array = self.removeDarkImage(image_array)
            return image_array * np.mean(self.flat_field_image_array) \
                   / (self.flat_field_image_array - self.dark_image_array)

    def removeDarkImage(self, image_array):
        """Uses dark image to correct fiber image

        Args:
            image_array: numpy array of the fiber image

        Returns:
            Corrected fiber image
        """
        image_array -= self.dark_image_array

        # Prevent any pixels from becoming negative values
        for x in xrange(self.image_width):
            for y in xrange(self.image_height):
                if image_array[y, x] < 0:
                    image_array[y, x] = 0.0

        return image_array

#==============================================================================#
#==== Private Variable Getters and Setters ====================================#
#==============================================================================#

    def getImageArray(self):
        return self.image_array

    def getImageHeight(self):
        return self.image_height

    def getImageWidth(self):
        return self.image_width

    def setImageArray(self, image):
        if type(image) is str:
            image = self.convertImageToArray(image)
        self.image_array = image

    def setDarkImage(self, *images):
        """Sets the corrective dark image

        Args:
            *images: either a single numpy array or multiple image file names

        Returns:
            None
        """
        if type(images[0]) is str:
            self.dark_image_array = np.zeros((self.image_height, self.image_width))
            for image_string in images:
                self.dark_image_array += self.convertImageToArray(image_string) \
                                         / float(len(images))

            # Remove any dark current from each pixel
            self.image_array -= self.dark_image_array

            for x in xrange(self.image_width):
                for y in xrange(self.image_height):
                    if self.image_array[y, x] < 0:
                        self.image_array[y, x] = 0.0

        else:
            self.dark_image_array = images[0]

    def setFlatFieldImage(self, *images):
        """Sets the corrective flat field image

        Args:
            *images: either a single numpy array or multiple image file names

        Returns:
            None
        """
        if type(images[0]) is str:
            self.flat_field_image_array = np.zeros((self.image_height, self.image_width))
            for image_string in images:
                self.flat_field_image_array += self.convertImageToArray(image_string) \
                                               / float(len(images))

        else:
            self.flat_field_image_array = images[0]

    def setBacklitImage(self, image):
        """Sets the backlit fiber image

        Args:
            *images: either a single numpy array or multiple image file names

        Returns:
            None
        """
        if type(image) is str:
            self.backlit_image_array = self.convertImageToArray(image_string)
        else:
            self.backlit_image_array = image

#==============================================================================#
#==== Image Centroiding =======================================================#
#==============================================================================#

    def getFiberCentroid(self):
        """Returns the centroid of the image

        Returns:
            (centroid_y, centroid_x)
        """
        if self.centroid_y is None or self.centroid_x is None:
            row_sum = self.getRowSum(self.image_array)
            column_sum = self.getColumnSum(self.image_array)

            row_weighted_sum = 0
            row_weight = 0
            for i in xrange(self.image_width):
                row_weighted_sum += i * row_sum[i]
                row_weight += row_sum[i]
            centroid_column = row_weighted_sum/row_weight

            columnWeightedSum = 0
            columnWeight = 0
            for i in xrange(self.image_height):
                columnWeightedSum += i * column_sum[i]
                columnWeight += column_sum[i]
            centroid_row = columnWeightedSum/columnWeight

            self.centroid_y = centroid_row
            self.centroid_x = centroid_column

        return self.centroid_y, self.centroid_x

    def getRowSum(self):
        row_sum = np.sum(self.image_array, axis=0)
        return ((row_sum - np.min(row_sum)) / self.image_height).astype(float)

    def getColumnSum(self):
        column_sum = np.sum(self.image_array, axis=1)
        return ((column_sum - np.min(column_sum)) / self.image_width).astype(float)

#==============================================================================#
#==== Image Centering =========================================================#
#==============================================================================#

    def getFiberCenterRadiusIteration(self, image_array=None):
        """Finds fiber center using a dark circle with various radii

        Uses a golden mean optimization method to find the optimal radius of the
        dark circle that covers the fiber image used in
        getFiberCenterCircleIteration(). The optimization is for a parameter
        array_sum which is weighted by the area of the circle, meaning that a
        smaller circle is preferred over one that simply covers the entire image

        Args:
            image_array: numpy array of the fiber image. If None, will use the
                         backlit_image_array

        Returns:
            (optimal_radius, center_y, center_x, optimal_array_sum)
        """
        if image_array is None:
            image_array = self.backlit_image_array

        r = np.zeros(4).astype(float)
        r[0] = 0
        r[3] = min(self.image_width, self.image_height) / 2
        r[1] = r[0] + (1 - self.phi) * (r[3] - r[0])
        r[2] = r[0] + self.phi * (r[3] - r[0])

        array_sum = np.zeros(2).astype(float)
        for i in xrange(2):
            center_y, center_x, array_sum[i] = \
                self.getFiberCenterCircleIteration(r[i+1], image_array)
            array_sum[i] *= r[i+1]

        while abs(r[3]-r[0]) > 1:
            min_index = np.argmin(array_sum) # Integer 0 or 1 for min of r[1], r[2]

            if min_index == 0:
                r[3] = r[2]
                r[2] = r[1]
                r[1] = r[0] + (1 - self.phi) * (r[3] - r[0])
            else:
                r[0] = r[1]
                r[1] = r[2]
                r[2] = r[0] + self.phi * (r[3] - r[0])

            array_sum[1 - min_index] = array_sum[min_index]

            center_y, center_x, array_sum[min_index] = \
                self.getFiberCenterCircleIteration(r[min_index+1], image_array)
            array_sum[min_index] *= r[min_index+1]

        return r[min_index+1], center_y, center_x, np.amin(array_sum)

    def getFiberCenterSquareIteration(self, radius, image_array=None):
        """Finds fiber center using a dark square of half-length radius

        Uses golden mean method to find the optimal center for a square covering
        the fiber image. The optimization is for a parameter array_sum that
        simply sums over the entire fiber image array

        Args:

        """

    def getFiberCenterCircleIteration(self, radius, image_array=None):
        """Finds fiber center using a dark circle of set radius

        Uses golden mean method to find the optimal center for a circle covering
        the fiber image. The optimization is for a parameter array_sum that
        simply sums over the entire fiber image array

        Args:
            radius: circle radius to test
            image_array: numpy array of the fiber image

        Returns:
            (center_y, center_x)
        """
        if image_array is None:
            image_array = self.backlit_image_array

        print "Testing Radius:", radius
        # Create four "corners" to test center of the removed circle
        x = np.zeros(4).astype(float)
        x[0] = radius
        x[3] = self.image_width - radius
        x[1] = x[0] + (1 - self.phi) * (x[3] - x[0])
        x[2] = x[0] + self.phi * (x[3] - x[0])

        y = np.zeros(4).astype(float)
        y[0] = radius
        y[3] = self.image_height - radius
        y[1] = y[0] + (1 - self.phi) * (y[3] - y[0])
        y[2] = y[0] + self.phi * (y[3] - y[0])

        # Initialize the array sum to each corner
        array_sum = np.zeros((2, 2)).astype(float)
        for i in xrange(2):
            for j in xrange(2):
                # Center at (y[j+1], x[i+1])
                removed_circle_array = self.removeCircle(radius, x[i+1], y[j+1], image_array)
                array_sum[j, i] = self.getArraySum(removed_circle_array)

        while abs(x[2] - x[1]) > 1 and abs(y[2] - y[1]) > 1:
            # Find the index of the minimum array sum corner
            min_index = np.unravel_index(np.argmin(array_sum), (2, 2)) # Tuple

            # Move the other corners to smaller search area
            if min_index[0] == 0:
                y[3] = y[2]
                y[2] = y[1]
                y[1] = y[0] + (1 - self.phi) * (y[3] - y[0])
            else:
                y[0] = y[1]
                y[1] = y[2]
                y[2] = y[0] + self.phi * (y[3] - y[0])
            if min_index[1] == 0:
                x[3] = x[2]
                x[2] = x[1]
                x[1] = x[0] + (1 - self.phi) * (x[3] - x[0])
            else:
                x[0] = x[1]
                x[1] = x[2]
                x[2] = x[0] + self.phi * (x[3] - x[0])

            # Replace the opposite corner array sum
            array_sum[1 - min_index[0], 1 - min_index[1]] = array_sum[min_index]

            # Recalculate new sums for other three corners
            for i in xrange(2):
                for j in xrange(2):
                    if i != 1 - min_index[1] or j != 1 - min_index[0]:
                        removed_circle_array = self.removeCircle(radius, x[i+1], y[j+1], image_array)
                        array_sum[j, i] = self.getArraySum(removed_circle_array)

        self.showImageArray(removed_circle_array)

        center_x = x[min_index[1]+1]
        center_y = y[min_index[0]+1]
        minarray_sum = np.amin(array_sum)

        return center_y, center_x, minarray_sum

    def getArraySum(self, image_array=None):
        if image_array is None:
            image_array = self.image_array
        return np.sum(image_array)

    def removeCircle(self, radius, x, y, image_array=None):
        """Removes a circle from an array

        Args:
            radius: circle's radius
            x: horizontal pixel number for circle center
            y: vertical pixel number for circle center
            image_array: numpy array of fiber image

        Returns:
            output_array: a copy of image_array with the circle removed
        """
        if image_array is None:
            image_array = self.image_array

        output_array = np.copy(image_array)
        circle_array = self.circleArray(radius, x, y, 1)
        output_array -= output_array * circle_array

        return output_array

    def circleArray(self, radius, x, y, res=1):
        """
            circle_array(radius, x, y)

         Returns:
            2D numpy array of integer values where points inside the circle are
            1 and outside the circle are 0. Points along the edge of the circle
            are weighted based on their relative distance to the center
        """
        x_int = int(round(x))
        y_int = int(round(y))
        radius_int = int(round(radius))

        circle_array = np.zeros((self.image_height, self.image_width)).astype(float)
        for i in xrange(x_int - radius_int - 1, x_int + radius_int + 2):
            for j in xrange(y_int - radius_int - 1, y_int + radius_int + 2):
                if (x-i)**2 + (y-j)**2 < (radius-1)**2:
                    circle_array[j, i] = 1.0
                elif (x-i)**2 + (y-j)**2 <= (radius+1)**2:
                    for k in xrange(res):
                        for l in xrange(res):
                            d = (x-(i - 0.5 + 1/(2.0*res) + k/(1.0*res)))**2 \
                                + (y-(j - 0.5 + 1/(2.0*res) + l/(1.0*res)))**2
                            if d <= radius**2:
                                circle_array[j, i] += (1.0/res)**2

        return circle_array

    def getFiberCenterEdgeMethod(self):
        """The averages of the fiber edges gives the fiber center

        Returns:
            center_y, center_x
        """
        if self.center_y_edge is None or self.center_x_edge is None:
            if self.left is None or self.top is None:
                self.setFiberEdges()

            self.center_y_edge = (self.top + self.bottom)/2.0
            self.center_x_edge = (self.left + self.right)/2.0

        return self.center_y_edge, self.center_x_edge

    def getFiberDiameter(self):
        if self.fiber_diameter is None:
            self.setFiberEdges()
        return self.fiber_diameter

    def setFiberEdges(self, threshold=None):
        """Set fiber edge pixel values

        Sets the left, right, top, and bottom edges of the fiber by finding where
        the row and column sums cross the given threshold. Also sets the width
        of the fiber by the maximum of the horizontal and vertical lengths

        Args:
            threshold: the minimum brightness value used to determine the fiber
                edge

        Returns:
            None
        """
        if threshold is None:
            threshold = self.threshold

        row_sum = self.getRowSum(self.image_array)
        column_sum = self.getColumnSum(self.image_array)

        left = -1
        right = -1
        for index, value in enumerate(row_sum):
            if left < 0:
                if value > threshold:
                    left = index
            else:
                if value > threshold:
                    right = index
        top = -1
        bottom = -1
        for index, value in enumerate(column_sum):
            if top < 0:
                if value > threshold:
                    top = index
            else:
                if value > threshold:
                    bottom = index

        self.left = left
        self.right = right
        self.top = top
        self.bottom = bottom
        self.fiber_diameter = max(right - left, bottom - top)

#=============================================================================#
#===== Image Plotting ========================================================#
#=============================================================================#

    def plotHorizontalCrossSection(self, row):
        plt.plot(self.image_array[row, :])
        plt.title('Horizontal Cross Section (row = %s)'%row)
        plt.xlabel('Pixel')
        plt.show()

    def plotVerticalCrossSection(self, column):
        plt.plot(self.image_array[:, column])
        plt.title('Vertical Cross Section (column = %s)'%column)
        plt.xlabel('Pixel')

    def plotCrossSections(self, row, column):
        plt.figure(1)
        plt.subplot(211)
        self.plotHorizontalCrossSection(row)
        plt.subplot(212)
        self.plotVerticalCrossSection(column)
        plt.show()

    def plotCrossSectionSums(self):
        plt.figure(1)
        plt.subplot(211)
        plt.plot(self.getRowSum())
        plt.title('Average for each Column')
        plt.xlabel('Column')
        plt.subplot(212)
        plt.plot(self.getColumnSum())
        plt.title('Average for each Row')
        plt.xlabel('Row')
        plt.show()

    def showImageArray(self, image_array=None):
        if image_array is None:
            image_array = self.image_array
        plt.figure(1)
        plt.imshow(image_array)
        plt.show()

    def iterationCounter(self, iteration):
        iteration += 1
        if iteration % 10 == 0:
            print iteration
        return iteration
