import numpy as np
import matplotlib.pyplot as plt
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

    def getImageData(self, image_obj):
        y0, x0 = image_obj.getFiberCenter(method='edge', tol=1, test_range=10, show_image=False)
        radius = image_obj.getFiberRadius()        
        image_array = image_obj.getImageArray()
        #image_array = image_obj.getTophatFit()
        return image_array, y0, x0, radius

#=============================================================================#
#==== Modal Noise Methods ====================================================#
#=============================================================================#

    def getModalNoise(self, image_obj=None, method='fft', output='array', radius_factor=0.95, deg=8):
        """Finds modal noise of camera image using specified method

        Args:
            camera: string designating near, far, or input
            method: string designating tophat, polynomial, gaussian, gradient,
                contrast, fft, gini, or entropy

        Returns:
            modal noise parameter
        """
        if image_obj is None:
            image_obj = self.image_obj

        if len(method) < 3:
            raise ValueError('Incorrect string for method type')
        elif method in 'tophat':
            return self.getModalNoiseTophat(image_obj, radius_factor)
        elif method in 'polynomial':
            return self.getModalNoisePolynomial(image_obj, radius_factor, deg=deg)
        elif method in 'gaussian':
            return self.getModalNoiseGaussian(image_obj, radius_factor)
        elif method in 'gradient':
            return self.getModalNoiseGradient(image_obj, radius_factor)
        elif method in 'contrast':
            return self.getModalNoiseContrast(image_obj, radius_factor)
        elif method in 'fft' or method in 'fourier':
            return self.getModalNoiseFFT(image_obj, radius_factor)
        elif method in 'gini':
            return self.getModalNoiseGini(image_obj, radius_factor)
        elif method in 'entropy':
            return self.getModalNoiseEntropy(image_obj, radius_factor)
        else:
            raise ValueError('Incorrect string for method type')

    def getModalNoiseTophat(self, image_obj=None, output='parameter', radius_factor=0.95):
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

        image_array, y0, x0 = self.cropImage(iamge_array, x0, y0, radius)
        
        if len(output) < 3:
            raise ValueError('Incorrect output string')

        elif output in 'array':
            return image_array

        elif output in 'parameter':
            intensity_array = self.getIntensityArray(image_array, x0, y0, radius)
            return intensity_array.std() / intensity_array.mean()

        else:
            raise ValueError('Incorrect output string')

    def getModalNoiseGradient(self, image_obj=None, output='array', radius_factor=0.95):
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

        if len(output) < 3:
            raise ValueError('Incorrect output string')

        elif output in 'array':
            return gradient_array

        elif output in 'parameter':
            intensity_array = self.getIntensityArray(gradient_array, x0, y0, radius)
            image_intensity_array = self.getIntensityArray(image_array, x0, y0, radius)
            return intensity_array.std() / image_intensity_array.mean()

        else:
            ValueError('Incorrect output string')

    def getModalNoiseFFT(self, image_obj=None, radius_factor=0.95):
        """Finds modal noise of image using the image's power spectrum
        
        Args:
            image_obj (optional): ImageAnalysis object to be used
            radius_factor (optional): fraction of the radius inside which the
                modal noise calculation will be made

        Returns:
            fft_list, freq_list  
        """
        if image_obj is None:
            image_obj = self.image_obj

        image_array, y0, x0, radius = self.getImageData(image_obj)
        height, width = image_array.shape

        image_array, x0, y0 = self.cropImage(image_array, x0, y0, min(x0, y0, width-x0, height-y0))
        height, width = image_array.shape
        image_array = self.isolateCircle(image_array, x0, y0, radius*1.0)

        #fft_array = np.abs(np.fft.fft2(image_array, norm='ortho')[:height/2, :width/2])
        fft_array = np.fft.fftshift(np.abs(np.fft.fft2(image_array, norm='ortho')))
        fx0 = width/2
        fy0 = height/2

        max_freq = 100 #min(height, width)/2
        bin_width = 0.05
        list_len = int(max_freq / bin_width) + 1

        fft_list = np.zeros(list_len).astype('float64')
        weight_list = np.zeros(list_len).astype('float64')
        freq_list = bin_width * np.arange(list_len).astype('float64')

        for i in range(fx0-max_freq, fx0+max_freq):
            for j in range(fy0-max_freq, fy0+max_freq):
                freq = np.sqrt((fx0-i)**2 + (fy0-j)**2)
                if freq <= max_freq:
                    fft_list[int(freq/bin_width)] += fft_array[j,i]
                    weight_list[int(freq/bin_width)] += 1.0

        mask = (weight_list > 0.0).astype('bool')
        weight_list = weight_list[mask]
        freq_list = freq_list[mask]
        fft_list = fft_list[mask] / weight_list # Average out

        fft_list /= fft_list.max() # Normalize
        freq_list /= width * image_obj.pixel_size / 1000 # Get frequencies in 1/mm

        if len(output) < 3:
            raise ValueError('Incorrect output string')

        elif output in 'array':
            return fft_list, freq_list

        elif output in 'parameter':
            return self.getGiniCoefficient(fft_list)

        else:
            ValueError('Incorrect output string')

    def getModalNoiseGini(self, image_obj=None, radius_factor=0.95):
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

        if len(output) < 3:
            raise ValueError('Incorrect output string')

        elif output in 'array':
            return self.cropImage(image_array, x0, y0, radius)[0]

        elif output in 'parameter':
            return self.getGiniCoefficient(intensity_array)
        else:
            raise ValueError('Incorrect output string')

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

    def getModalNoiseContrast(self, image_obj=None, radius_factor=0.95):
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

        if len(output) < 3:
            raise ValueError('Incorrect output string')

        elif output in 'array':
            return    

        return (intensity_array.max() - intensity_array.min()) / (intensity_array.max() + intensity_array.min())

    def getModalNoisePolynomial(self, image_obj=None, radius_factor=0.95, deg=8):
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

    def getModalNoiseGaussian(self, image_obj=None, radius_factor=0.95):
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

    def getModalNoiseEntropy(self, image_obj=None, radius_factor=0.95):
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

if __name__ == '__main__':
    from ImageAnalysis import ImageAnalysis
    import os as os
    import matplotlib.pyplot as plt
    plt.rc('font', size=16, family='sans-serif')
    plt.rc('xtick', labelsize=14)
    plt.rc('ytick', labelsize=14)
    plt.rc('lines', lw=2)

    calibration_folder = "Calibration/TIF/"
    image_folder = '../Alignment Images/2016-07-12/'

    nf_led_images = []
    for i in xrange(3):
        nf_led_images.append(image_folder + 'nf_led_0um_' + str(i) + '_80ms.tif')

    nf_laser_images = []
    nf_laser_images_agitated = []
    for i in xrange(10):
        nf_laser_images.append(image_folder + 'nf_laser_noagitation_' + str(i) + '_1.8ms.tif')
        nf_laser_images_agitated.append(image_folder + 'nf_laser_agitation_' + str(i) + '_1.8ms.tif')    

    nf_dark_images = []
    nf_ambient_images = []
    for i in xrange(3):
        nf_dark_images.append(calibration_folder + 'nf_dark_' + str(i) + '_10ms.tif')
        nf_ambient_images.append(image_folder + 'nf_ambient_' + str(i) + '_1.8ms.tif')

    nf_flat_images = []
    #for i in xrange(8):
    #    nf_flat_images.append(calibration_folder + 'nf_flat_' + str(i) + '_1ms.tif')

    LED_nf = ImageAnalysis(nf_led_images, nf_dark_images, nf_flat_images, nf_ambient_images, 3.45, 10, 16)
    print 'LED image initialized'
    laser_nf_agitated = ImageAnalysis(nf_laser_images_agitated, nf_dark_images, nf_flat_images, nf_ambient_images, 3.45, 10, 16)
    print 'Agitated laser image initialized'
    laser_nf = ImageAnalysis(nf_laser_images, nf_dark_images, nf_flat_images, nf_ambient_images, 3.45, 10, 16)
    print 'Unagitated laser image initialized'
    print

    laser_mn = ModalNoise(laser_nf)
    laser_mn_agitated = ModalNoise(laser_nf_agitated)
    LED_mn = ModalNoise(LED_nf)

    LED_fft, LED_freq = LED_mn.getModalNoise(method='fft')
    print 'LED finished'
    laser_ag_fft, laser_ag_freq = laser_mn_agitated.getModalNoise(method='fft')
    print 'Agitated laser finished'
    laser_fft, laser_freq = laser_mn.getModalNoise(method='fft')
    print 'Unagitated laser finished'

    plt.figure(1)
    plt.plot(LED_freq, LED_fft, label='LED')
    plt.plot(laser_ag_freq, laser_ag_fft, label='Agitated laser')
    plt.plot(laser_freq, laser_fft, label='Unagitated laser')
    plt.xlim(0, 10)
    plt.yscale('log')
    plt.ylabel('Normalized Power')
    plt.xlabel('Frequency [1/mm]')
    plt.legend()
    plt.show()
