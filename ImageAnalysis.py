"""ImageAnalysis.py was written by Ryan Petersburg for use with fiber
characterization on the EXtreme PRecision Spectrograph
"""
import numpy as np
from NumpyArrayHandler import *
from Calibration import Calibration

class ImageAnalysis(object):
    """Fiber face image analysis class

    Class that conducts image analysis on a fiber face image after it has been
    corrected by the given dark and flat field images. Also contains information
    about the CCD that took the image. Public methods in this class allow
    calculation of the image's centroid as well as multiple methods to find the
    fiber center and diameter

    Attributes:
        image
        pixel_size
        camera
        kernel_size
        magnification
        test
        exp_time
        bit_depth
        date_time
        temp
        num_images
    """
    def __init__(self, image_input, calibration=None,
                 pixel_size=None, camera=None, threshold=1000,
                 kernel_size=3, magnification=None):
    	self._image_info = dict(height=None,
    		 				    width=None,
    						    pixel_size=pixel_size,
    						    camera=camera,
    						    magnification=magnification,
    						    exp_time=None,
    						    bit_depth=None,
    						    date_time=None,
    						    temp=None,
    						    num_images=None)

    	self._analysis_info = dict(kernel_size=kernel_size,
    		                       threshold=threshold)

        self._calibration = calibration
        if self._calibration is None:
            self._calibration = Calibration(None, None, None)

        self.image = None
        self.setImageArray(image_input)
        self._filtered_image = self.getFilteredImage(self.image, self.kernel_size)

        self.threshold = threshold
        if self.threshold is None:
            mask = (self._filtered_image > self._calibration.dark_image.max()).astype('bool')
            self.threshold = self._filtered_image[mask].mean() / 2.0
        
        self._edges = dict(left=None,
                           right=None,
                           top=None,
                           bottom=None)
        self._center = dict(edge=dict(x=None, y=None),
                            radius=dict(x=None, y=None),
                            circle=dict(x=None, y=None),
                            gaussian=dict(x=None, y=None))
        self._centroid = dict(x=None, y=None)
        self._diameter = dict(edge=None,
                              radius=None
                              circle=None
                              gaussian=None)
        self._array_sum = dict(radius=None,
                               circle=None)

        self._fit = dict(gaussian=None,
        				 polynomial=None)

        # Golden Ratio for optimization tests
        self._phi = (5 ** 0.5 - 1) / 2

#=============================================================================#
#==== Private Variable Setters ===============================================#
#=============================================================================#

    def setImageArray(self, image_input):
        """Sets the uncorrected image to be analyzed

        Args:
            image_input: see convertImageToArray for options

        Sets:
            self.image
        """
        self._uncorrected_image, output_dict = convertImageToArray(image_input, True)

        self.setImageInfo(output_dict)

        self.image = self._calibration.executeErrorCorrections(self._uncorrected_image,
                                                                self.image_info['exp_time'])

    def setImageInfo(self, output_dict):
        self.image_info['height'], self.image_info['width'] = self.image.shape

    	for key in output_dict:
    		self.image_info[key] = output_dict[key]

        if 'pixel_size' in output_dict:
            self.image_info['pixel_size'] = output_dict['pixel_size']
        if 'bit_depth' in output_dict:
            self.image_info['bit_depth'] = output_dict['bit_depth']
        if 'exp_time' in output_dict:
            self.image_info['exp_time'] = output_dict['exp_time']
        if 'date_time' in output_dict:
            self.image_info['date_time'] = output_dict['date_time']
        if 'temp' in output_dict:
            self.image_info['temp'] = output_dict['temp']
        if 'camera' in output_dict:
            self.image_info['camera'] = output_dict['camera']
        if 'num_images' in output_dict:
            self.image_info['num_images'] = output_dict['num_images']

        if self.magnification is None:
            if self.image_info['camera'] == 'nf' or self.image_info['camera'] == 'in':
                self.magnification = 10.0
            if self.image_info['camera'] == 'ff':
                self.magnification = 1.0

        if self.image_info['pixel_size'] is None:
            raise RuntimeError('Pixel Size needs to be set externally')
        if self.magnification is None:
            raise RuntimeError('Magnification needs to be set externally')

    def saveData(self, folder=None):
        data = dict(image_info=self.image_info,
                    image=list(self.image),
                    )
        
        with open(folder + 'ImageAnalysisData.txt', 'w') as file:
            file.write

#=============================================================================#
#==== Private Variable Getters ===============================================#
#=============================================================================#

    def getHeight(self):
        """Getter for the image height

        Returns:
            self.image_info['height']
        """
        return self.image_info['height']

    def getWidth(self):
        """Getter for the image width

        Returns:
            self.image_info['width']
        """
        return self.image_info['width']

    def getFiberData(self, method=None, show_image=False, tol=1, test_range=None, units='microns'):
        """Getter for the fiber center and diameter

        Returns:
            (center_y, center_x, diameter)
        """
        return (self.getFiberCenter(method, show_image, tol, test_range, units) 
                + (self.getFiberDiameter(method, show_image, tol, test_range, units),))

    def getFiberRadius(self, method=None, show_image=False, tol=1, test_range=None, units='pixels'):
        """Getter for the fiber radius

        Finds the radius of the fiber using the given method or, if no method
        is given, the most precise method already completed

        Args:
            method (optional): string representing the fiber centering method

        Returns:
            fiber radius
        """
        return self.getFiberDiameter(method, show_image, tol, test_range, units) / 2.0

    def getFiberDiameter(self, method=None, show_image=False, tol=1,
                         test_range=None, units='pixels'):
        """Getter for the fiber diameter in pixels

        Find the diameter of the fiber using the given method or, if no method
        is given, the most precise method already completed

        Args:
            method (optional): string representing the fiber centering method

        Returns:
            fiber diameter
        """
        if method is None:
            if self._diameter['radius'] is not None:
                method = 'radius'
            elif self._diameter['gaussian'] is not None:
                method = 'gaussian'
            else:
                method = 'edge'

        if self._diameter[method] is None:
        	self.setFiberDiameter(method, tol, test_range)

        diameter = self._diameter[method]

        if units == 'pixels':
            return diameter
        elif units == 'microns':
            return diameter * self.image_info['pixel_size'] / self.magnification
        else:
            raise RuntimeError('Incorrect string for units')

    def getFiberCenter(self, method=None, show_image=False, tol=1, test_range=None, units='pixels'):
        """Getter for the fiber center in pixels

        Find the center position of the fiber using the given method or, if no
        method is given, the most precise method already completed

        Args:
            method (optional): string representing the fiber centering method
            show_image (optional): boolean for whether or not to show image of
                completed method
            tol (optional): tolerance value passed to
                getFiberCenterRadiusMethod()and getFiberCenterCircleMethod()
            test_range (optional): range of tested values passed to
                getFiberCenterRadiusMethod() and getFiberCenterCircleMethod()

        Returns:
            center y, center x
        """
        if method is None:
            if self._center['radius']['x'] is not None:
                method = 'radius'
            elif self._center['gaussian']['x'] is not None:
                method = 'gaussian'
            elif self._center['circle']['x'] is not None:
                method = 'circle'
            else:
                method = 'edge'

        if self._center[method]['x'] is None:
        	self.setFiberCenter(method, tol, test_range)

        if show_image:
        	self.showOverlaidTophat(self._center[method]['x'],
        		                    self._center[method]['y'],
        		                    self._diameter[method] / 2.0,
        		                    tol=1)

        center = self._center[method]['y'], self._center[method]['x']

        if units == 'pixels':
        	return center
        elif units == 'microns':
        	return tuple(np.array(center) * self.image_info['pixel_size'] / self.magnification)
        else:
        	raise RuntimeError('Incorrect string for units')

    def getFiberCenterRadiusMethod(self, tol=1, test_range=None, show_image=False):
        """Getter for the fiber center using the radius method

        See setFiberCenterRadiusMethod() for method details

        Args:
            show_image (optional): boolean for whether or not to show image of
                completed method

        Returns:
            center y (pixels), center x (pixels)
        """
        if self._center['radius']['x'] is None:
            self.setFiberCenterRadiusMethod(tol, test_range)
            
        if show_image:
            self.showOverlaidTophat(self._center['radius']['x'],
                                    self._center['radius']['y'],
                                    self._diameter['radius'] / 2.0,
                                    tol=1)

        return self._center['radius']['y'], self._center['radius']['x']

    def getFiberCenterCircleMethod(self, radius=None, tol=1, test_range=None, show_image=False):
        """Getter for the fiber center using the circle method

        See setFiberCenterCircleMethod() for method details

        Args:
            show_image (optional): boolean for whether or not to show image of
                completed method

        Returns:
            center y (pixels), center x (pixels)
        """
        if radius is None:
            radius = self.getFiberRadius()
        if self._center['circle']['x'] is not None:
            show_image = False

        self.setFiberCenterCircleMethod(radius, tol, test_range)

        if show_image:
            self.showOverlaidTophat(self._center['circle']['x'],
                                    self._center['circle']['y'],
                                    self._diameter['circle'] / 2.0,
                                    tol=1)

        return self._center['circle']['y'], self._center['circle']['x']

    def getFiberCenterEdgeMethod(self, show_image=False):
        """Getter for the fiber center using the edge method

        See setFiberCenterEdgeMethod() for method details

        Args:
            show_image (optional): boolean for whether or not to show image of
                completed method

        Returns:
            center y (pixels), center x (pixels)
        """
        if self._center['edge']['y'] is None:
            self.setFiberCenterEdgeMethod()

        if show_image:
            self.showOverlaidTophat(self._center['edge']['x'],
                                    self._center['edge']['y'],
                                    self._diameter['edge'] / 2.0)

        return self._center['edge']['y'], self._center['edge']['x']

    def getFiberCenterGaussianMethod(self, show_image=False):
        """Getter for the fiber center using the gaussian method

        See setFiberCenterGaussianMethod() for method details

        Args:
            show_image (optional): boolean for whether or not to show image of
                completed method

        Returns:
            center y (pixels), center x (pixels)
        """
        if self._center['gaussian']['x'] is None:
            self.setFiberCenterGaussianMethod()

        if show_image:
            showImageArray(self._fit['gaussian'])
            plotOverlaidCrossSections(self.image, self._fit['gaussian'],
                                      self._center['gaussian']['y'], self._center['gaussian']['x'])

        return self._center['gaussian']['y'], self._center['gaussian']['x']

    def getFiberCentroid(self, radius_factor=None, method='edge', units='pixels'):
        """Getter for the fiber centroid

        See setFiberCentroid() for method details

        Returns:
            centroid y (pixels), centroid x (pixels)
        """
        if self._centroid['x'] is None:
            self.setFiberCentroid(radius_factor, method)
        centroid = (self._centroid['y'], self._centroid['x'])

        if units == 'microns':
            return tuple(np.array(centroid) * self.image_info['pixel_size'] / self.magnification)
        return centroid

    def getGaussianFit(self, image_array=None, initial_guess=None, full_output=False):
        if image_array is None:
            if self._fit['gaussian'] is None:
                self.setFiberCenterGaussianMethod()
            return self._fit['gaussian']
        return gaussianFit(image_array, initial_guess, full_output)

    def getMeshGrid(self, image_array=None):
        if image_array is None:
            image_array = self.image
        return meshGridFromArray(image_array)

    def getPolynomialFit(self, image_array=None, deg=6, x0=None, y0=None):
        if image_array is None:
            image_array = self.image
        return polynomialFit(image_array, deg, x0, y0)

    def getTophatFit(self):
        y0, x0 = self.getFiberCenter(show_image=False)
        radius = self.getFiberRadius(show_image=False)
        return circleArray(self.getMeshGrid(), x0, y0, radius, res=1)

    def getFilteredImage(self, image_array=None, kernel_size=None):
        if image_array is None and kernel_size is None and self._filtered_image is not None:
            return self._filtered_image
        if image_array is None:
            image_array = self.image
        if kernel_size is None:
            kernel_size = self.kernel_size
        return filteredImage(image_array, kernel_size)

    def getDarkImage(self):
        return calibration.dark_image

    def getAmbientImage(self):
        return calibration.ambient_image

    def getFlatImage(self):
        return calibration.flat_image

    def getArraySum(self, image_array=None):
        if image_array is None:
            image_array = self.image
        sumArray(image_array)

#=============================================================================#
#==== Image Centroiding ======================================================#
#=============================================================================#

    def setFiberCentroid(self, radius_factor=None, method='edge'):
        """Finds the centroid of the fiber face image

        Args:
            radius_factor: the factor by which the radius is multiplied when
                isolating the fiber face in the image

        Sets:
            centroid_y
            centroid_x
        """
        if radius_factor is not None:
            y0, x0 = self.getFiberCenter(method=method, show_image=False)
            radius = self.getFiberRadius(method=method, show_image=False)
            image_array_iso = isolateCircle(self.image, x0, y0,
                                            radius*radius_factor, res=1)
        else:
            image_array_iso = self.image

        x_array, y_array = self.getMeshGrid()
        self._centroid['x'] = (image_array_iso * x_array).sum() / image_array_iso.sum()
        self._centroid['y'] = (image_array_iso * y_array).sum() / image_array_iso.sum()

#=============================================================================#
#==== Image Centering ========================================================#
#=============================================================================#

	def setFiberDiameter(self, method, tol=1, test_range=None):
		"""Finds fiber diameter using given method

		See each respective method for centering algorithm

		Raises:
			RuntimeError: cannot accept the 'circle' method when setting the
				diameter since it requires a known radius to run
		"""
		if method == 'circle':
			raise RuntimeError('Fiber diameter cannot be set by circle method')
		self.setFiberCenter(method, tol, test_range)

	def setFiberCenter(self, method, tol=1, test_range=None, radius=None):
		"""Find fiber center using given method
		
		See each respective method for centering algorithm

		Raises:
			RuntimeError: needs a valid method string to run the proper
				algorithm
		"""
		if method == 'radius':
			self.setFiberCenterRadiusMethod(tol, test_range)
		elif method == 'edge':
			self.setFiberCenterEdgeMethod()
		elif method == 'circle':
			self.setFiberCenterCircleMethod(radius, tol, test_range)
		elif method == 'gaussian':
			self.setFiberCenterGaussianMethod()
		else:
			raise RuntimeError('Incorrect string for fiber centering method')

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
        y0, x0 = self.getFiberCenter(method='edge', show_image=False)
        initial_guess = (x0, y0, self.getFiberRadius(),
                         self.image.max(), self.image.min())

        self._fit['gaussian'], opt_parameters = self.getGaussianFit(self._filtered_image,
                                                                 initial_guess=initial_guess,
                                                                 full_output=True)

        self._center['gaussian']['x'] = opt_parameters[0]
        self._center['gaussian']['y'] = opt_parameters[1]
        self._diameter['gaussian'] = opt_parameters[2] * 2

    def setFiberCenterRadiusMethod(self, tol=1, test_range=None):
        """Finds fiber center using a dark circle with various radii

        Uses a golden mean optimization method to find the optimal radius of the
        dark circle that covers the fiber image used in
        getFiberCenterCircleMethod(). The optimization is for a parameter
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
            if r[0] < 0.0:
                r[0] = 0.0
            r[3] = approx_radius + test_range
        else:
            r[0] = 0
            r[3] = min(self.image_info['height'], self.image_info['width']) / 2.0

        r[1] = r[0] + (1 - self._phi) * (r[3] - r[0])
        r[2] = r[0] + self._phi * (r[3] - r[0])

        array_sum = np.zeros(2).astype(float)
        for i in xrange(2):
            self.setFiberCenterCircleMethod(r[i+1], tol, test_range)
            array_sum[i] = self._array_sum['circle'] + self.threshold * np.pi * r[i+1]**2

        min_index = np.argmin(array_sum) # Integer 0 or 1 for min of r[1], r[2]

        while abs(r[3]-r[0]) > tol:
            if min_index == 0:
                r[3] = r[2]
                r[2] = r[1]
                r[1] = r[0] + (1 - self._phi) * (r[3] - r[0])
            else:
                r[0] = r[1]
                r[1] = r[2]
                r[2] = r[0] + self._phi * (r[3] - r[0])

            array_sum[1 - min_index] = array_sum[min_index]

            self.setFiberCenterCircleMethod(r[min_index+1], tol, test_range)
            array_sum[min_index] = (self._array_sum['circle']
                                    + self.threshold * np.pi * r[min_index+1]**2)

            min_index = np.argmin(array_sum) # Integer 0 or 1 for min of r[1], r[2]

        self._diameter['radius'] = r[min_index+1] * 2
        self._center['radius']['y'] = self._center['circle']['y']
        self._center['radius']['x'] = self._center['circle']['x']
        self._array_sum['radius'] = np.amin(array_sum)

    def setFiberCenterCircleMethod(self, radius, tol=1, test_range=None):
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
        res = 1 #int(1.0/tol)

        # Create four "corners" to test center of the removed circle
        x = np.zeros(4).astype(float)
        y = np.zeros(4).astype(float)

        if test_range is not None:
            approx_center = self.getFiberCenter(method='edge', show_image=False)
            test_range = test_range / 2.0

            x[0] = approx_center[1] - test_range
            if x[0] < radius:
                x[0] = radius
            x[3] = approx_center[1] + test_range
            if x[3] > self.image_info['width'] - radius:
                x[3] = self.image_info['width'] - radius

            y[0] = approx_center[0] - test_range
            if y[0] < radius:
                y[0] = radius
            y[3] = approx_center[0] + test_range
            if y[3] > self.image_info['height'] - radius:
                y[3] = self.image_info['height'] - radius

        else:
            x[0] = radius
            x[3] = self.image_info['width'] - radius

            y[0] = radius
            y[3] = self.image_info['height'] - radius

        x[1] = x[0] + (1 - self._phi) * (x[3] - x[0])
        x[2] = x[0] + self._phi * (x[3] - x[0])

        y[1] = y[0] + (1 - self._phi) * (y[3] - y[0])
        y[2] = y[0] + self._phi * (y[3] - y[0])

        # Initialize array sums to each corner
        array_sum = np.zeros((2, 2)).astype(float)
        for i in xrange(2):
            for j in xrange(2):
                removed_circle_array = removeCircle(self._filtered_image,
                                                    x[i+1], y[j+1],
                                                    radius, res)
                array_sum[j, i] = sumArray(removed_circle_array)

        # Find the index of the corner with minimum array_sum
        min_index = np.unravel_index(np.argmin(array_sum), (2, 2)) # Tuple

        while abs(x[3] - x[0]) > tol and abs(y[3] - y[0]) > tol:
            # Move the other corners to smaller search area
            if min_index[0] == 0:
                y[3] = y[2]
                y[2] = y[1]
                y[1] = y[0] + (1 - self._phi) * (y[3] - y[0])
            else:
                y[0] = y[1]
                y[1] = y[2]
                y[2] = y[0] + self._phi * (y[3] - y[0])
            if min_index[1] == 0:
                x[3] = x[2]
                x[2] = x[1]
                x[1] = x[0] + (1 - self._phi) * (x[3] - x[0])
            else:
                x[0] = x[1]
                x[1] = x[2]
                x[2] = x[0] + self._phi * (x[3] - x[0])

            # Replace the opposite corner array sum (so it doesn't need to be recalculated)
            array_sum[1 - min_index[0], 1 - min_index[1]] = array_sum[min_index]
            min_index = (1 - min_index[0], 1 - min_index[1])

            # Recalculate new sums for all four corners
            for i in xrange(2):
                for j in xrange(2):
                    if i != min_index[1] or j != min_index[0]:
                        removed_circle_array = removeCircle(self._filtered_image,
                                                            x[i+1], y[j+1],
                                                            radius, res)
                        array_sum[j, i] = sumArray(removed_circle_array)

            min_index = np.unravel_index(np.argmin(array_sum), (2, 2))

        self._center['circle']['x'] = x[min_index[1]+1]
        self._center['circle']['y'] = y[min_index[0]+1]
        self._diameter['circle'] = radius * 2.0
        self._array_sum['circle'] = np.amin(array_sum)

    def setFiberCenterEdgeMethod(self):
        """The averages of the fiber edges gives the fiber center

        Returns:
            center_y, center_x
        """
        self.setFiberEdges()

        self._center['edge']['y'] = (self._edges['top'] + self._edges['bottom']) / 2.0
        self._center['edge']['x'] = (self._edges['left'] + self._edges['right']) / 2.0

    def setFiberEdges(self):
        """Set fiber edge pixel values

        Sets the left, right, top, and bottom edges of the fiber by finding where
        the maxima of each row and column cross the given threshold. Also sets
        the width of the fiber by the maximum of the horizontal and vertical
        lengths

        Sets:
            self._edges['left']
            self._edges['right']
            self._edges['top']
            self._edges['bottom']
            self._diameter['edge']
        """
        left = -1
        right = -1
        for index in xrange(self.image_info['width']):
            if left < 0:
                if self._filtered_image[:, index].max() > self.threshold:
                    left = index
            else:
                if self._filtered_image[:, index].max() > self.threshold:
                    right = index

        top = -1
        bottom = -1
        for index in xrange(self.image_info['height']):
            if top < 0:
                if self._filtered_image[index, :].max() > self.threshold:
                    top = index
            else:
                if self._filtered_image[index, :].max() > self.threshold:
                    bottom = index

        self._edges['left'] = left
        self._edges['right'] = right
        self._edges['top'] = top
        self._edges['bottom'] = bottom
        self._diameter['edge'] = ((right - left) + (bottom - top)) / 2.0

#=============================================================================#
#==== Overriding Methods =====================================================#
#=============================================================================#
    def plotCrossSections(self, image_array=None, row=None, column=None):
        if image_array is None:
            image_array = self.image
        if row is None:
            row = self.image_info['height'] / 2.0
        if column is None:
            column = self.image_info['width'] / 2.0
        plotCrossSections(image_array, row, column)

    def showImageArray(self, image_array=None):
        if image_array is None:
            image_array = self.image
        showImageArray(image_array)


    def showOverlaidTophat(self, x0, y0, radius, tol=1):
        res = int(1.0/tol)
        showImageArray(removeCircle(self.image, x0, y0, radius, res=res))
        plotOverlaidCrossSections(self._filtered_image,
                                  2 * self.threshold
                                  *circleArray(self.getMeshGrid(),
                                               x0, y0, radius, res=res),
                                  y0, x0)


if __name__ == "__main__":
    folder = 'Scrambling Measurements/Core Extension/2016-08-05 Prototype Core Extension 1/'

    calibration = Calibration([folder + 'Dark/in_' + str(i).zfill(3) + '.fit' for i in xrange(10)],
                              None,
                              [folder + 'Ambient/in_' + str(i).zfill(3) + '.fit' for i in xrange(10)])

    nf_images = [folder + 'Shift_00/in_' + str(i).zfill(3) + '.fit' for i in xrange(10)]

    im_obj = ImageAnalysis(nf_images, calibration, threshold=200, magnification=10, kernel_size=11)

    tol = 1
    test_range = 50
    factor = 1.0

    print 'Height:', im_obj.getHeight(), 'Width:', im_obj.getWidth()
    #im_obj.showImageArray()
    print 'Kernel:', im_obj.kernel_size
    print 'Threshold:', im_obj.threshold
    print
    print 'Centroid'
    centroid_row, centroid_column = im_obj.getFiberCentroid(factor)
    print 'Centroid Row:', centroid_row, 'Centroid Column:', centroid_column
    print
    print 'Edge:'
    center_y, center_x = im_obj.getFiberCenter(method='edge', show_image=True)
    print 'Diameter:', im_obj.getFiberDiameter(method='edge', show_image=False, units='microns'), 'microns'
    print 'Center Row:', center_y, 'Center Column:', center_x
    print
    print 'Radius:'
    center_y, center_x = im_obj.getFiberCenter(method= 'radius', tol=tol, show_image=True, test_range=test_range)
    print 'Diameter:', im_obj.getFiberDiameter(method='radius',show_image=False, units='microns'), 'microns'
    print 'Center Row:', center_y, 'Center Column:', center_x
    print
    print 'Gaussian:'
    center_y, center_x = im_obj.getFiberCenter(method='gaussian')
    print 'Diameter:', im_obj.getFiberDiameter(method='gaussian', units='microns'), 'microns'
    print 'Center Row:', center_y, 'Center Column:', center_x
