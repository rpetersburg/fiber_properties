import numpy as np
from scipy import optimize as opt
from NumpyArrayHandler import NumpyArrayHandler

class ImageAnalysis(NumpyArrayHandler):

    def __init__(self, image_string, dark_image_strings, flat_image_strings, pixel_size=5.2, magnification=10):
        self.threshold = 1
        self.pixel_size = pixel_size / magnification
        self.bit_rate = 8

        self.setImageArray(image_string)
        self.image_height, self.image_width = self.image_array.shape

        self.setDarkImage(*dark_image_strings)
        self.setFlatFieldImage(*flat_image_strings)

        self.executeErrorCorrections()

        #Approximate the Fiber Center and Diameter
        self.setFiberCenterEdgeMethod()

        self.setMeshGrid()

        self.centroid_y = None
        self.centroid_x = None

        self.center_y_circle = None
        self.center_x_circle = None
        self.center_y_radius = None
        self.center_x_radius = None
        self.center_y_gaussian = None
        self.center_x_gaussian = None
        self.center_y_tophat = None
        self.center_x_tophat = None

        self.fiber_diameter_radius = None
        self.fiber_diameter_gaussian = None
        self.fiber_diameter_tophat = None

        self.gaussian_fit = None
        self.tophat_fit = None

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

        elif self.fiber_diameter_tophat is not None:
            return self.getFiberDiameterTophatMethod()

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

    def getFiberDiameterTophatMethod(self):
        if self.fiber_diameter_tophat is None:
            self.setFiberCenterTophatMethod()            
        return self.fiber_diameter_tophat

    def getFiberCenter(self):
        if self.center_x_gaussian is not None:
            return self.getFiberCenterGaussianMethod()

        elif self.center_x_radius is not None:
            return self.getFiberCenterRadiusIteration()

        elif self.center_x_circle is not None:
            return self.getFiberCenterCircleIteration()

        return self.center_y_edge, self.center_x_edge

    def getFiberCenterRadiusIteration(self):
        if self.center_x_radius is None:
            self.setFiberCenterRadiusIteration()
        return self.center_y_radius, self.center_x_radius

    def getFiberCenterCircleIteration(self, radius=None):
        if radius is None:
            radius = self.getFiberRadius()
        if self.center_x_circle is None:
            self.setFiberCenterCircleIteration(radius)
        return self.center_y_circle, self.center_x_circle

    def getFiberCenterEdgeMethod(self):
        return self.center_y_edge, self.center_x_edge

    def getFiberCenterGaussianMethod(self):
        if self.center_x_gaussian is None:
            self.setFiberCenterGaussianMethod()
        return self.center_y_gaussian, self.center_x_gaussian

    def getFiberCenterTophatMethod(self):
        if self.center_x_tophat is None:
            self.setFiberCenterTophatMethod()
        return self.center_y_tophat, self.center_x_tophat

    def getFiberCentroid(self):
        if self.centroid_x is None:
            self.setFiberCentroid()
        return self.centroid_y, self.centroid_x

    def getGaussianFit(self):
        if self.gaussian_fit is None:
            self.setGaussianFit()

        return self.gaussian_fit

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
        self.mesh_grid = np.meshgrid(np.linspace(0, self.image_width - 1,self.image_width),
                                     np.linspace(0, self.image_height - 1, self.image_height))

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

        opt_parameters, cov_matrix = opt.curve_fit(self.gaussianArray, (self.mesh_grid),
                                                   self.image_array.ravel(), p0=initial_guess,
                                                   method='dogbox')

        self.gaussian_fit = self.gaussianArray(self.mesh_grid, *opt_parameters).reshape(self.image_height,
                                                                                          self.image_width)

        self.center_x_gaussian = opt_parameters[0]
        self.center_y_gaussian = opt_parameters[1]
        self.fiber_diameter_gaussian = opt_parameters[2] * 2

        self.showImageArray(self.gaussian_fit)
        self.plotOverlaidCrossSections(self.image_array, self.gaussian_fit,
                                       self.center_y_gaussian, self.center_x_gaussian)

    def setFiberCenterTophatMethod(self):
        """Finds fiber center using a Tophat Fit
        
        Sets:
            fiber_diameter_tophat
            center_y_tophat
            center_x_tophat
        """
        approx_center_x = self.getFiberCenter()[1]
        approx_center_y = self.getFiberCenter()[0]
        approx_radius = self.getFiberRadius()
        approx_amp = self.image_array.max()

        initial_guess = (approx_center_x, approx_center_y, approx_radius)

        bounds = [(approx_center_x - 10, approx_center_x + 10),
                  (approx_center_y - 10, approx_center_y + 10),
                  (approx_radius - 10, approx_radius + 10)]

        opt_parameters = opt.differential_evolution(self.removedCircleArraySum, bounds=bounds, tol=1).x

        self.tophat_fit = self.circleArray((self.mesh_grid), *opt_parameters) * 256

        self.center_x_tophat = opt_parameters[0]
        self.center_y_tophat = opt_parameters[1]
        self.fiber_diameter_tophat = opt_parameters[2] * 2

        #self.showImageArray(self.tophat_fit)
        self.plotOverlaidCrossSections(self.image_array, self.tophat_fit,
                                       self.center_y_tophat, self.center_x_tophat)

    def removedCircleArraySum(self, var):
        """
        Args:
            var[0]: x0
            var[1]: y0
            var[2]: radius
        """
        return self.getArraySum(self.removeCircle(var[2], var[0], var[1], self.image_array)) * var[2]**2

    def setFiberCenterRadiusIteration(self):
        """Finds fiber center using a dark circle with various radii

        Uses a golden mean optimization method to find the optimal radius of the
        dark circle that covers the fiber image used in
        getFiberCenterCircleIteration(). The optimization is for a parameter
        array_sum which is weighted by the area of the circle, meaning that a
        smaller circle is preferred over one that simply covers the entire image

        Sets:
            fiber_diameter_radius: diameter of the fiber (radius method)
            center_y_radius: y-position of center (radius method)
            center_x_radius: x-position of center (radius method)
        """
        r = np.zeros(4).astype(float)
        r[0] = 0
        r[3] = min(self.image_width, self.image_height) / 2
        r[1] = r[0] + (1 - self.phi) * (r[3] - r[0])
        r[2] = r[0] + self.phi * (r[3] - r[0])

        array_sum = np.zeros(2).astype(float)
        for i in xrange(2):
            self.setFiberCenterCircleIteration(r[i+1])
            array_sum[i] = self.array_sum_circle * r[i+1]**2

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

            self.setFiberCenterCircleIteration(r[min_index+1])
            array_sum[min_index] = self.array_sum_circle * r[min_index+1]**2

        self.fiber_diameter_radius = r[min_index+1] * 2
        self.center_y_radius = self.center_y_circle
        self.center_x_radius = self.center_x_circle
        self.array_sum_radius = np.amin(array_sum)

    def setFiberCenterCircleIteration(self, radius):
        """Finds fiber center using a dark circle of set radius

        Uses golden mean method to find the optimal center for a circle
        covering the fiber image. The optimization is for a parameter array_sum
        that simply sums over the entire fiber image array

        Args:
            radius: circle radius to test
        """
        #print "Testing Radius:", radius

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
                removed_circle_array = self.removeCircle(radius, x[i+1], y[j+1], self.image_array)
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
                        removed_circle_array = self.removeCircle(radius, x[i+1], y[j+1], self.image_array)
                        array_sum[j, i] = self.getArraySum(removed_circle_array)

        #self.showImageArray(removed_circle_array)

        self.center_x_circle = x[min_index[1]+1]
        self.center_y_circle = y[min_index[0]+1]
        self.array_sum_circle = np.amin(array_sum)

    def removeCircle(self, radius, x0, y0, image_array):
        """Removes a circle from an array

        Args:
            radius: circle's radius
            x: horizontal pixel number for circle center
            y: vertical pixel number for circle center

        Returns:
            output_array: a copy of image_array with the circle removed
        """
        output_array = np.copy(image_array)
        circle_array = self.circleArray(self.mesh_grid, x0, y0, radius)
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
        self.fiber_diameter_edge = max(right - left, bottom - top)

#=============================================================================#
#==== Overriding Methods =====================================================#
#=============================================================================#

    def showImageArray(self, image_array=None):
        if image_array is None:
            image_array = self.image_array
        super(ImageAnalysis, self).showImageArray(image_array)