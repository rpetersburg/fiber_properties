from NumpyArrayHandler import NumpyArrayHandler


class ModalNoise(NumpyArrayHandler):

    def __init__(self, image_obj, camera='nf'):
        self.setImageObject(image_obj)
        self.setCamera(camera)        

    def setImageObject(self, image_obj):
        self.image_obj = image_obj

    def setCamera(self, camera):
        self.camera = camera

    def getImageObject(self):
        return self.image_obj

    def getCamera(self):
        return self.camera

#=============================================================================#
#==== Modal Noise Methods ====================================================#
#=============================================================================#

    def getModalNoise(self, image_obj=None, method='fft', radius_factor=0.9, deg=8):
        """Finds modal noise of camera image using specified method

        Args:
            camera: string designating near, far, or input
            method: string designating tophat, polynomial, or gaussian

        Returns:
            modal noise parameter
        """
        if image_obj is None:
            image_obj = self.image_obj

        if method == 'tophat':
            return self.getModalNoiseTophat(image_obj, radius_factor)
        elif method == 'polynomial':
            return self.getModalNoisePolynomial(image_obj, radius_factor, deg=deg)
        elif method == 'gaussian':
            return self.getModalNoiseGaussian(image_obj, radius_factor)
        elif method == 'gradient':
            return self.getModalNoiseGradient(image_obj, radius_factor)
        elif method == 'contrast':
            return self.getModalNoiseContrast(image_obj, radius_factor)
        elif method == 'fft':
            return self.getModalNoiseFFT(image_obj, radius_factor)
        elif method == 'gini':
            return self.getModalNoiseGini(image_obj, radius_factor)
        elif method == 'entropy':
            return self.getModalNoiseEntropy(image_obj, radius_factor)
        else:
            raise ValueError('Incorrect string for method type')

    def getModalNoiseTophat(self, image_obj=None, radius_factor=0.9):
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
            image_obj = self.image_obj

        image_array, y0, x0, radius = self.getImageData(image_obj)
        radius *= radius_factor

        intensity_array = self.getIntensityArray(image_array, x0, y0, radius)

        image_array, x0, y0 = self.cropImage(image_array, x0, y0, radius)
        self.plotCrossSections(image_array, y0, x0)

        return intensity_array.std() / intensity_array.mean()

    def getModalNoiseGradient(self, image_obj=None, radius_factor=0.9):
        """Finds modal noise of image using the image gradient

        Args:
            image_obj (optional): ImageAnalysis object to be used
            radius_factor (optional): fraction of the radius inside which the
                modal noise calculation will be made

        Returns:
            modal noise parameter
        """
        if image_obj is None:
            image_obj = self.image_obj

        image_array, y0, x0, radius = self.getImageData(image_obj)
        image_array, x0, y0 = self.cropImage(image_array, x0, y0, radius)
        radius *= radius_factor

        gradient_y, gradient_x = np.gradient(image_array)
        gradient_array = np.sqrt(gradient_x**2 + gradient_y**2)

        self.showImageArray(gradient_array)
        #self.plotOverlaidCrossSections(image_array, gradient_array, y0, x0)

        intensity_array = self.getIntensityArray(gradient_array, x0, y0, radius)
        image_intensity_array = self.getIntensityArray(image_array, x0, y0, radius)

        return intensity_array.std() / image_intensity_array.mean()

    def getModalNoiseFFT(self, image_obj=None, radius_factor=0.9):
        """Finds modal noise of image using the image's power spectrum
        
        Args:
            image_obj (optional): ImageAnalysis object to be used
            radius_factor (optional): fraction of the radius inside which the
                modal noise calculation will be made        
        """
        if image_obj is None:
            image_obj = self.image_obj

        image_array, y0, x0, radius = self.getImageData(image_obj)
        height, width = image_array.shape

        image_array, x0, y0 = self.cropImage(image_array, width/2, height/2, min(height, width)/2 - 1)
        height, width = image_array.shape
        self.showImageArray(image_array)

        fft_array = np.abs(np.fft.fft2(image_array, norm='ortho')[:height/2, :width/2])

        f_list = []
        fft_list = []
        for i in xrange(min(height, width)/2):
            for j in xrange(i+1):
                if np.sqrt(i**2 + j**2) < 70:
                    fft_list.append((fft_array[j,i] + fft_array[i,j]) / 2.0)
                    f_list.append(np.sqrt(i**2 + j**2))

        f_list = np.array(f_list).astype('float64')*1000/3.45/width
        fft_list = np.array(fft_list).astype('float64') / sum(fft_list)
        print f_list
        print fft_list

        plt.figure(1)
        plt.scatter(f_list, fft_list)
        plt.show()

        #self.showImageArray(fft_array)
        #self.plotCrossSections(fft_array, 0, 0)

        return #1 - self.getGiniCoefficient(fft_array)

    def getModalNoiseGini(self, image_obj=None, radius_factor=0.9):
        """Find modal noise of image using Gini coefficient

        Args:
            image_obj (optional): ImageAnalysis object to be used
            radius_factor (optional): fraction of the radius inside which the
                modal noise calculation will be made

        Returns:
            modal noise parameter
        """
        if image_obj is None:
            image_obj = self.image_obj

        image_array, y0, x0, radius = self.getImageData(image_obj)
        radius *= radius_factor

        intensity_array = self.getIntensityArray(image_array, x0, y0, radius)
        return self.getGiniCoefficient(intensity_array)

    def getGiniCoefficient(self, test_array):
        """Finds gini coefficient for intensities in given array

        Args:
            test_array: arbitrary-dimension numpy array

        Returns:
            gini coefficient
        """
        test_array = test_array.ravel()
        gini_coeff = 0.0
        for i in xrange(test_array.size):
            gini_coeff += np.abs(test_array[i] - test_array).sum()
        norm = 2.0 * test_array.sum() * test_array.size

        return gini_coeff / norm

    def getModalNoiseEntropy(self, image_obj=None, radius_factor=0.9):
        """Find modal noise of image using hartley entropy

        Args:
            image_obj (optional): ImageAnalysis object to be used
            radius_factor (optional): fraction of the radius inside which the
                modal noise calculation will be made

        Returns:
            modal noise parameter
        """
        if image_obj is None:
            image_obj = self.image_obj

        image_array, y0, x0, radius = self.getImageData(image_obj)
        radius *= radius_factor

        intensity_array = self.getIntensityArray(image_array, x0, y0, radius)
        intensity_array = intensity_array / intensity_array.sum()
        return (-intensity_array * np.log10(intensity_array)).sum()

    def getModalNoiseContrast(self, image_obj=None, radius_factor=0.9):
        """Finds modal noise of image using Michelson contrast
        
        Modal noise is defined as (I_max - I_min) / (I_max + I_min)

        Args:
            image_obj (optional): ImageAnalysis object to be used
            radius_factor (optional): fraction of the radius inside which the
                modal noise calculation will be made

        Returns:
            modal noise parameter
        """
        if image_obj is None:
            image_obj = self.image_obj

        image_array, y0, x0, radius = self.getImageData(image_obj)
        radius *= radius_factor

        image_array, x0, y0 = self.cropImage(image_array, x0, y0, radius)
        self.plotCrossSections(image_array, y0, x0)

        intensity_array = self.getIntensityArray(image_array, x0, y0, radius)        

        return (intensity_array.max() - intensity_array.min()) / (intensity_array.max() + intensity_array.min())

    def getModalNoisePolynomial(self, image_obj=None, radius_factor=0.9, deg=8):
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
            image_obj = self.image_obj

        image_array, y0, x0, radius = self.getImageData(image_obj)
        radius *= radius_factor

        image_array, x0, y0 = self.cropImage(image_array, x0, y0, radius)

        poly_fit = self.getPolynomialFit(image_array, deg=deg, x0=x0, y0=y0)
        #self.showImageArray(poly_fit)
        self.plotOverlaidCrossSections(image_array, poly_fit, radius, radius)

        diff_array = image_array - poly_fit

        intensity_array = self.getIntensityArray(diff_array, x0, y0, radius)
        image_intensity_array = self.getIntensityArray(image_array, x0, y0, radius)

        return intensity_array.std() / image_intensity_array.mean()

    def getModalNoiseGaussian(self, image_obj=None, radius_factor=0.9):
        """Finds modal noise of image using gaussian fit
        
        Modal noise is defined as the variance in the difference between the
        original image and the gaussian fit normalized by the mean of the
        intensity inside the fiber face

        Args:
            image_obj (optional): ImageAnalysis object to be used
            radius_factor (optional): fraction of the radius inside which the
                modal noise calculation will be made

        Returns:
            modal noise parameter
        """
        if image_obj is None:
            image_obj = self.image_obj

        image_array, y0, x0, radius = self.getImageData(image_obj)
        radius *= radius_factor

        image_array, x0, y0 = self.cropImage(image_array, x0, y0, radius)
        gaussian_fit = self.getGaussianFit(image_array)

        diff_array = image_array - gaussian_fit
        #self.showImageArray(gaussian_fit)
        self.plotOverlaidCrossSections(image_array, gaussian_fit, y0, x0)

        intensity_array = self.getIntensityArray(diff_array, x0, y0, radius)
        image_intensity_array = self.getIntensityArray(image_array, x0, y0, radius)

        return intensity_array.std() / image_intensity_array.mean()

    def getImageData(self, image_obj):
        y0, x0 = image_obj.getFiberCenter(method='edge', tol=1, test_range=10, show_image=True)
        radius = image_obj.getFiberRadius()        
        image_array = image_obj.getImageArray()
        #image_array = image_obj.getGaussianFit()
        return image_array, y0, x0, radius