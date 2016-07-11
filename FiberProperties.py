import numpy as np
from NumpyArrayHandler import NumpyArrayHandler

class FiberProperties(NumpyArrayHandler):

    def __init__(self, in_object=None, nf_object=None, ff_object=None):
        self.in_object = in_object
        self.nf_object = nf_object
        self.ff_object = ff_object
        self.nf_scrambling_gain = None
        self.ff_scrambling_gain = None

#=============================================================================#
#==== Private Variable Setters ===============================================#
#=============================================================================#

    def setInputObject(self, in_object):
        self.in_object = in_object

    def setNearFieldObject(self, nf_object):
        self.nf_object = nf_object

    def setFarFieldObject(self, ff_object):
        self.ff_object = ff_object

#=============================================================================#
#==== Scrambling Gain Methods ================================================#
#=============================================================================#

    def getNearFieldScramblingGain(self):
        if self.nf_scrambling_gain is None:
            self.nf_scrambling_gain = getScramblingGain(self.in_object, self.nf_object)
        return self.nf_scrambling_gain

    def getFarFieldScramblingGain(self):
        if self.ff_scrambling_gain is None:
            self.ff_scrambling_gain = getScramblingGain(self.in_object, self.ff_object)
        return self.ff_scrambling_gain

    def getScramblingGain(self, in_object, out_object):
        in_centroid_y, in_centroid_x = in_object.getFiberCentroid()
        in_y0, in_x0 = in_object.getFiberCenterEdgeMethod()
        in_diameter = in_object.getFiberDiameter()

        out_centroid_y, out_centroid_x = out_object.getFiberCentroid()
        out_y0, out_x0 = out_object.getFiberCenterEdgeMethod()
        out_diameter = out_object.getFiberDiameter()

        delta_D_in = np.sqrt((in_centroid_x - in_x0)**2 + (in_centroid_y - in_y0)**2)
        delta_D_out = np.sqrt((out_centroid_x - out_x0)**2 + (out_centroid_y - out_y0)**2)

        scramblingGain = (delta_D_in / in_diameter) / (delta_D_out / out_diameter)

        return scramblingGain

#=============================================================================#
#==== Modal Noise Methods ====================================================#
#=============================================================================#

    def getModalNoise(self, camera, method, radius_factor=0.8):
        """Finds modal noise of camera image using specified method

        Args:
            camera: string designating near, far, or input
            method: string designating tophat, polynomial, or gaussian

        Returns:
            modal noise parameter
        """
        if camera == 'near' or camera == 'nf':
            image_obj = self.nf_object
        elif camera == 'far' or camera == 'ff':
            image_obj = self.ff_object
        elif camera == 'input' or camera == 'in':
            image_obj = self.in_object
        else:
            raise ValueError('Incorrect string for camera type')

        if method == 'tophat':
            return self.getModalNoiseTophat(image_obj, radius_factor)
        elif method == 'polynomial':
            return self.getModalNoisePolynomial(image_obj, radius_factor)
        elif method == 'gaussian':
            return self.getModalNoiseGaussian(image_obj, radius_factor)
        elif method == 'gradient':
            return self.getModalNoiseGradient(image_obj, radius_factor)
        else:
            raise ValueError('Incorrect string for method type')

    def getModalNoiseTophat(self, image_obj=None, radius_factor=0.95):
        """Finds modal noise of image assumed to be a tophat
        
        Modal noise is defined as the variance across the fiber face normalized
        by the mean across the fiber face

        Args:
            image_obj (optional): ImageAnalysis object to be used
            radius_factor (optional): fraction of the radius inside which the
                modal noise calculation will be made

        Returns:
            modal noise parameter
        """
        if image_obj is None:
            image_obj = self.nf_object

        image_array = image_obj.getImageArray()
        y0, x0 = image_obj.getFiberCenterRadiusIteration(tol=1, test_range=30)
        radius = image_obj.getFiberRadius() * radius_factor

        intensity_array = self.getIntensityArray(image_array, x0, y0, radius)

        return intensity_array.std() / intensity_array.mean()

    def getModalNoiseGradient(self, image_obj=None, radius_factor=0.95):
        """Finds modal noise of image using the image gradient

        Args:
            image_obj (optional): ImageAnalysis object to be used
            radius_factor (optional): fraction of the radius inside which the
                modal noise calculation will be made

        Returns:
            modal noise parameter
        """
        if image_obj is None:
            image_obj = self.nf_object

        image_array = image_obj.getImageArray()
        y0, x0 = image_obj.getFiberCenterRadiusIteration(tol=1, test_range=30)
        radius = image_obj.getFiberRadius() * radius_factor

        gradient_y, gradient_x = np.gradient(image_array)
        gradient_array = np.sign(gradient_x) * np.sqrt(gradient_x**2 + gradient_y**2)

        self.showImageArray(gradient_array)
        self.plotOverlaidCrossSections(image_array, gradient_array, y0, x0)

        intensity_array = self.getIntensityArray(gradient_array, x0, y0, radius)
        image_intensity_array = self.getIntensityArray(image_array, x0, y0, radius)

        return intensity_array.std() / image_intensity_array.mean()

    def getModalNoisePolynomial(self, image_obj=None, radius_factor=0.95):
        """Finds modal noise of image using polynomial fit
        
        Crops image exactly around the circumference of the circle and fits a
        polynomial to the 2D image. Modal noise is then defined as the variance
        in the difference between the original image and the polynomial fit
        normalized by the maximum value of the polynomial fit

        Args:
            image_obj (optional): ImageAnalysis object to be used
            radius_factor (optional): fraction of the radius inside which the
                modal noise calculation will be made

        Returns:
            modal noise parameter
        """
        if image_obj is None:
            image_obj = self.ff_object

        image_array = image_obj.getImageArray()
        y0, x0 = image_obj.getFiberCenterRadiusIteration(tol=1, test_range=30)
        radius = image_obj.getFiberRadius() * radius_factor

        image_crop = image_array[int(y0-radius):int(y0+radius)+2,
                                 int(x0-radius):int(x0+radius)+2]
        y0 = radius + (y0-radius)-int(y0-radius)
        x0 = radius + (x0-radius)-int(x0-radius)
        
        poly_fit = self.getPolynomialFit(image_crop, deg=4, x0=x0, y0=y0)
        self.showImageArray(poly_fit)
        self.plotOverlaidCrossSections(image_crop, poly_fit, radius, radius)

        diff_array = image_crop - poly_fit
        self.plotCrossSections(diff_array, radius, radius)

        intensity_array = self.getIntensityArray(diff_array, x0, y0, radius)
        poly_intensity_array = self.getIntensityArray(poly_fit, x0, y0, radius)

        return intensity_array.std() / poly_intensity_array.mean()

    def getModalNoiseGaussian(self, image_obj=None, radius_factor=1.0):
        """Finds modal noise of image using gaussian fit
        
        Modal noise is defined as the variance in the difference between the
        original image and the gaussian fit normalized by the maximum value of
        the polynomial fit

        Args:
            image_obj (optional): ImageAnalysis object to be used
            radius_factor (optional): fraction of the radius inside which the
                modal noise calculation will be made

        Returns:
            modal noise parameter
        """
        if image_obj is None:
            image_obj = self.ff_object

        image_array = image_obj.getImageArray()
        gaussian_fit = image_obj.getGaussianFit()
        y0, x0 = image_obj.getFiberCenterGaussianMethod()
        radius = image_obj.getFiberRadius() * radius_factor

        diff_array = image_array - gaussian_fit
        self.plotCrossSections(diff_array, y0, x0)

        intensity_array = self.getIntensityArray(diff_array, x0, y0, radius)
        gaussian_intensity_array = self.getIntensityArray(gaussian_fit, x0, y0, radius)

        return intensity_array.std() / gaussian_intensity_array.mean()

    def getIntensityArray(self, image_array, x0, y0, radius):
        """Finds intensities inside a circle

        Returns an array of intensities from image_array which are contained 
        within the circle centered at (x0, y0) with radius radius

        Returns:
            intensity_array: one-dimensional numpy array of intensities
        """
        height, width = image_array.shape

        intensity_list = []
        for x in xrange(width):
            for y in xrange(height):
                if (x0-x)**2 + (y0-y)**2 <= (radius)**2:
                    intensity_list.append(image_array[y,x])

        return np.array(intensity_list)

#=============================================================================#
#==== Show Methods ===========================================================#
#=============================================================================#

    def showInputImageArray(self):
        self.showImageArray(self.in_object.getImageArray())

    def showNearFieldImageArray(self):
        self.showImageArray(self.nf_object.getImageArray())

    def showFarFieldImageArray(self):
        self.showImageArray(self.ff_object.getImageArray())

    def showImageArray(self, image_array=None):
        if image_array is not None:
            super(FiberProperties, self).showImageArray(image_array)
        else:
            if self.in_object is not None:
                self.showInputImageArray()
            if self.nf_object is not None:
                self.showNearFieldImageArray()
            if self.ff_object is not None:
                self.showFarFieldImageArray()
