import numpy as np
from scipy import optimize as opt
from NumpyArrayHandler import NumpyArrayHandler
import matplotlib.pyplot as plt

class ImageAnalysis(NumpyArrayHandler):

    def __init__(self, image_string, dark_image_strings, flat_image_strings, pixel_size=5.2, magnification=10):
        self.threshold = 1
        self.pixel_size = pixel_size
        self.magnification = magnification
        self.bit_rate = 8

        self.setImageArray(image_string)
        self.image_height, self.image_width = self.image_array.shape
        self.setMeshGrid()

        self.setDarkImage(*dark_image_strings)
        self.setFlatFieldImage(*flat_image_strings)

        self.executeErrorCorrections()

        #Approximate the Fiber Center and Diameter
        self.setFiberCenterEdgeMethod()        

        self.centroid_y = None
        self.centroid_x = None

        self.center_y_circle = None
        self.center_x_circle = None

        self.fiber_diameter_radius = None
        self.center_y_radius = None
        self.center_x_radius = None

        self.fiber_diameter_gaussian = None
        self.center_y_gaussian = None
        self.center_x_gaussian = None

        self.gaussian_fit = None

        # Golden Ratio for optimization tests
        self.phi = (5 ** 0.5 - 1) / 2

#=============================================================================#
#==== Image Initialization ===================================================#
#=============================================================================#

    def executeErrorCorrections(self):
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
        if self.flat_image_array is None or self.dark_image_array is None:
            if self.flat_image_array is None:
                print 'Flat Field Image has not been set'
            if self.dark_image_array is None:
                print 'Dark Image has not been set'
            return
        self.image_array = self.removeDarkImage(self.image_array)
        self.flat_image_array = self.removeDarkImage(self.flat_image_array)
        self.image_array = self.image_array * np.mean(self.flat_image_array) / self.flat_image_array

    def removeDarkImage(self, image_array):
        """Uses dark image to correct fiber image

        Args:
            image_array: numpy array of the fiber image

        Returns:
            image_array: corrected fiber image
        """
        image_array -= self.dark_image_array

        # Prevent any pixels from becoming negative values
        for x in xrange(self.image_width):
            for y in xrange(self.image_height):
                if image_array[y, x] <= 1.0:
                    image_array[y, x] = 0.0

        return image_array

#=============================================================================#
#==== Private Variable Getters ===============================================#
#=============================================================================#

    def getImageArray(self):
        return self.image_array

    def getImageHeight(self):
        return self.image_height

    def getImageWidth(self):
        return self.image_width

    def getFiberRadius(self):
        return self.getFiberDiameter() / 2

    def getFiberDiameter(self):
        if self.fiber_diameter_gaussian is not None:
            return self.getFiberDiameterGaussianMethod()

        elif self.fiber_diameter_radius is not None:
            return self.getFiberDiameterRadiusMethod()

        return self.getFiberDiameterEdgeMethod()

    def getFiberDiameterEdgeMethod(self):
        return self.fiber_diameter_edge

    def getFiberDiameterRadiusMethod(self):
        return self.fiber_diameter_radius

    def getFiberDiameterGaussianMethod(self):
        if self.fiber_diameter_gaussian is None:
            self.setFiberCenterGaussianMethod()
        return self.fiber_diameter_gaussian

    def getFiberCenter(self):
        if self.center_x_gaussian is not None:
            return self.getFiberCenterGaussianMethod()

        elif self.center_x_radius is not None:
            return self.getFiberCenterRadiusIteration()

        elif self.center_x_circle is not None:
            return self.getFiberCenterCircleIteration()

        return self.center_y_edge, self.center_x_edge

    def getFiberCenterRadiusIteration(self, tol=1, test_range=None):
        if self.center_x_radius is None:
            self.setFiberCenterRadiusIteration(tol, test_range)

        self.showOverlaidTophat(self.center_x_radius, self.center_y_radius, self.fiber_diameter_radius / 2.0)

        return self.center_y_radius, self.center_x_radius

    def getFiberCenterCircleIteration(self, radius=None, tol=1, test_range=None):
        if radius is None:
            radius = self.getFiberRadius()
        if self.center_x_circle is None:
            self.setFiberCenterCircleIteration(radius, tol, test_range)

        #self.showOverlaidTophat(self.center_x_circle, self.center_y_circle, self.fiber_diameter_circle / 2.0)

        return self.center_y_circle, self.center_x_circle

    def getFiberCenterEdgeMethod(self):
        self.showOverlaidTophat(self.center_x_edge, self.center_y_edge, self.fiber_diameter_edge / 2.0)

        return self.center_y_edge, self.center_x_edge

    def getFiberCenterGaussianMethod(self):
        if self.center_x_gaussian is None:
            self.setFiberCenterGaussianMethod()

        self.showImageArray(self.gaussian_fit)
        self.plotOverlaidCrossSections(self.image_array, self.gaussian_fit,
                                       self.center_y_gaussian, self.center_x_gaussian)

        return self.center_y_gaussian, self.center_x_gaussian

    def getFiberCentroid(self):
        if self.centroid_x is None:
            self.setFiberCentroid()
        return self.centroid_y, self.centroid_x

    def getGaussianFit(self, image_array=None, initial_guess=None, full_output=False):
        if image_array is None:
            if self.gaussian_fit is None:
                self.setFiberCenterGaussianMethod()
            return self.gaussian_fit
        return super(ImageAnalysis, self).getGaussianFit(image_array, initial_guess, full_output)

    def getMeshGrid(self, image_array=None):
        if image_array is None:
            return self.mesh_grid
        return super(ImageAnalysis, self).getMeshGrid(image_array)

    def getPolynomialFit(self, image_array=None, deg=3):
        if image_array is None:
            image_array = self.image_array
        return super(ImageAnalysis, self).getPolynomialFit(image_array, deg)

#=============================================================================#
#==== Private Variable Setters ===============================================#
#=============================================================================#

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
            self.dark_image_array = np.zeros((self.image_height,
                                              self.image_width))
            for image_string in images:
                self.dark_image_array += self.convertImageToArray(image_string) / float(len(images))
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
            self.flat_image_array = np.zeros((self.image_height,
                                                    self.image_width))
            for image_string in images:
                self.flat_image_array += self.convertImageToArray(image_string) / float(len(images))

        else:
            self.flat_image_array = images[0]

    def setBacklitImage(self, *images):
        """Sets the backlit fiber image

        Args:
            *images: either a single numpy array or multiple image file names

        Returns:
            None
        """
        if type(images[0]) is str:
            self.backlit_image_array = np.zeros((self.image_height, self.image_width))

            for image_string in images:
                self.backlit_image_array += self.convertImageToArray(image_string) / float(len(images))
        else:
            self.backlit_image_array = image

    def setMeshGrid(self):
        """Sets a two dimensional mesh grid

        Creates mesh_grid for the x and y pixels using the image height and width

        Sets:
            mesh_grid
        """
        self.mesh_grid = self.getMeshGrid(self.image_array)

#=============================================================================#
#==== Image Centroiding ======================================================#
#=============================================================================#

    def setFiberCentroid(self):
        """Finds the centroid of the image

        Sets:
            centroid_y
            centroid_x
        """
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

#=============================================================================#
#==== Image Centering ========================================================#
#=============================================================================#

    def setFiberCenterGaussianMethod(self):
        """Finds fiber center using a Gaussian Fit

        Uses Scipy.optimize.curve_fit method to fit fiber image to 
        self.gaussianArray. The radius found extends to 2-sigma of the gaussian
        therefore encompassing ~95% of the imaged light. Use previous methods
        of center-finding to approximate the location of the center

        Sets:
            fiber_diameter_gaussian: diameter of the fiber (gaussian method)
            center_y_gaussian: y-position of center (gaussian method)
            center_x_gaussian: x-position of center (gaussian method)
        """
        #initial_guess = (50,50,50,50)
        initial_guess = (self.getFiberCenter()[1], self.getFiberCenter()[0],
                         self.getFiberRadius(), self.image_array.max())

        self.gaussian_fit, opt_parameters = self.getGaussianFit(self.image_array, initial_guess, True)

        self.center_x_gaussian = opt_parameters[0]
        self.center_y_gaussian = opt_parameters[1]
        self.fiber_diameter_gaussian = opt_parameters[2] * 2

    def setFiberCenterRadiusIteration(self, tol, test_range):
        """Finds fiber center using a dark circle with various radii

        Uses a golden mean optimization method to find the optimal radius of the
        dark circle that covers the fiber image used in
        getFiberCenterCircleIteration(). The optimization is for a parameter
        array_sum which is weighted by the area of the circle, meaning that a
        smaller circle is preferred over one that simply covers the entire image

        Args:
            tol: minimum possible range of radius values before ending iteration
            test_range: range of tested radii. If None, uses full possible range

        Sets:
            fiber_diameter_radius: diameter of the fiber (radius method)
            center_y_radius: y-position of center (radius method)
            center_x_radius: x-position of center (radius method)
        """
        # Initialize range of tested radii
        r = np.zeros(4).astype(float)

        if test_range is not None:
            approx_radius = self.getFiberRadius()
            test_range = test_range / 2.0
        
            r[0] = approx_radius - test_range
            if r[0] < 0.0: r[0] = 0.0
            r[3] = approx_radius + test_range
        else:
            r[0] = 0
            r[3] = min(self.image_height, self.image_width) / 2.0

        r[1] = r[0] + (1 - self.phi) * (r[3] - r[0])
        r[2] = r[0] + self.phi * (r[3] - r[0])

        array_sum = np.zeros(2).astype(float)
        for i in xrange(2):
            self.setFiberCenterCircleIteration(r[i+1], tol, test_range)
            array_sum[i] = self.array_sum_circle + self.threshold * np.pi * r[i+1]**2

        min_index = np.argmin(array_sum) # Integer 0 or 1 for min of r[1], r[2]

        while abs(r[3]-r[0]) > tol:
            if min_index == 0:
                r[3] = r[2]
                r[2] = r[1]
                r[1] = r[0] + (1 - self.phi) * (r[3] - r[0])
            else:
                r[0] = r[1]
                r[1] = r[2]
                r[2] = r[0] + self.phi * (r[3] - r[0])

            array_sum[1 - min_index] = array_sum[min_index]

            self.setFiberCenterCircleIteration(r[min_index+1], tol, test_range)
            array_sum[min_index] = self.array_sum_circle + self.threshold * np.pi * r[min_index+1]**2

            min_index = np.argmin(array_sum) # Integer 0 or 1 for min of r[1], r[2]

        self.fiber_diameter_radius = r[min_index+1] * 2
        self.center_y_radius = self.center_y_circle
        self.center_x_radius = self.center_x_circle
        self.array_sum_radius = np.amin(array_sum)

    def setFiberCenterCircleIteration(self, radius, tol, test_range):
        """Finds fiber center using a dark circle of set radius

        Uses golden mean method to find the optimal center for a circle
        covering the fiber image. The optimization is for a parameter array_sum
        that simply sums over the entire fiber image array

        Args:
            radius: circle radius to test
            tol: minimum possible range of center_x or center_y values before
                ending iteration
            test_range: initial range of tested center values. If None, uses
                full range.
        """
        #print "Testing Radius:", radius
        res = int(1/tol)

        # Create four "corners" to test center of the removed circle
        x = np.zeros(4).astype(float)
        y = np.zeros(4).astype(float)

        if test_range is not None:
            approx_center = self.getFiberCenter()
            test_range = test_range / 2.0        
        
            x[0] = approx_center[1] - test_range
            if x[0] < radius: x[0] = radius
            x[3] = approx_center[1] + test_range
            if x[3] > self.image_width - radius:
                x[3] = self.image_width - radius

            y[0] = approx_center[0] - test_range
            if y[0] < radius: y[0] = radius
            y[3] = approx_center[0] + test_range
            if y[3] > self.image_height - radius:
                y[3] = self.image_height - radius
        
        else:
            x[0] = radius
            x[3] = self.image_width - radius

            y[0] = radius
            y[3] = self.image_height - radius

        x[1] = x[0] + (1 - self.phi) * (x[3] - x[0])
        x[2] = x[0] + self.phi * (x[3] - x[0])       
        
        y[1] = y[0] + (1 - self.phi) * (y[3] - y[0])
        y[2] = y[0] + self.phi * (y[3] - y[0])

        # Initialize array sums to each corner
        array_sum = np.zeros((2, 2)).astype(float)
        for i in xrange(2):
            for j in xrange(2):
                removed_circle_array = self.removeCircle(self.image_array, x[i+1], y[j+1], radius, res)
                array_sum[j, i] = self.getArraySum(removed_circle_array)

        # Find the index of the corner with minimum array_sum
        min_index = np.unravel_index(np.argmin(array_sum), (2, 2)) # Tuple

        while abs(x[3] - x[0]) > tol and abs(y[3] - y[0]) > tol:
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

            # Replace the opposite corner array sum (so it doesn't need to be recalculated)
            array_sum[1 - min_index[0], 1 - min_index[1]] = array_sum[min_index]
            min_index = (1 - min_index[0], 1 - min_index[1])

            # Recalculate new sums for all four corners
            for i in xrange(2):
                for j in xrange(2):
                    if i != min_index[1] or j != min_index[0]:
                        removed_circle_array = self.removeCircle(self.image_array, x[i+1], y[j+1], radius, res)
                        array_sum[j, i] = self.getArraySum(removed_circle_array)

            min_index = np.unravel_index(np.argmin(array_sum), (2, 2))

        self.center_x_circle = x[min_index[1]+1]
        self.center_y_circle = y[min_index[0]+1]
        self.fiber_diameter_circle = radius * 2.0
        self.array_sum_circle = np.amin(array_sum)

    def removeCircle(self, image_array, x0, y0, radius, res=1):
        """Removes a circle from an array

        Args:
            radius: circle's radius
            x: horizontal pixel number for circle center
            y: vertical pixel number for circle center

        Returns:
            output_array: a copy of image_array with the circle removed
        """
        output_array = np.copy(image_array)
        circle_array = self.circleArray(self.mesh_grid, x0, y0, radius, res)
        output_array -= output_array * circle_array

        return output_array

    def setFiberCenterEdgeMethod(self):
        """The averages of the fiber edges gives the fiber center

        Returns:
            center_y, center_x
        """
        self.setFiberEdges()

        self.center_y_edge = (self.top + self.bottom)/2.0
        self.center_x_edge = (self.left + self.right)/2.0

    def setFiberEdges(self):
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
        row_sum = self.getRowSum(self.image_array)
        column_sum = self.getColumnSum(self.image_array)

        left = -1
        right = -1
        for index, value in enumerate(row_sum):
            if left < 0:
                if value > self.threshold:
                    left = index
            else:
                if value > self.threshold:
                    right = index
        top = -1
        bottom = -1
        for index, value in enumerate(column_sum):
            if top < 0:
                if value > self.threshold:
                    top = index
            else:
                if value > self.threshold:
                    bottom = index

        self.left = left
        self.right = right
        self.top = top
        self.bottom = bottom
        self.fiber_diameter_edge = max(right - left, bottom - top)

#=============================================================================#
#==== Overriding Methods =====================================================#
#=============================================================================#

    def showImageArray(self, image_array=None):
        if image_array is None:
            image_array = self.image_array
        super(ImageAnalysis, self).showImageArray(image_array)


    def showOverlaidTophat(self, x0, y0, radius):
        self.showImageArray(self.removeCircle(self.image_array, x0, y0, radius))
        self.plotOverlaidCrossSections(self.image_array,
                                       self.image_array.max()*self.circleArray(self.mesh_grid, x0, y0, radius),
                                       y0, x0)