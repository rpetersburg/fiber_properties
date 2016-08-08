"""FiberProperties.py was written by Ryan Petersburg for use with fiber
characterization for the EXtreme PRecision Spectrograph

This method contains functions that calculate the scrambling gain, focal ratio
degradation, and modal noise for image analysis objects
"""

import numpy as np
import matplotlib.pyplot as plt

def getImageData(image_obj, method=None):
    """Returns relevant information from an ImageAnalysis object
    
    Args:
        image_obj [ImageAnalysis]: image object to be analyzed

    Returns:
        image_array [np.ndarray]: the 2D image
        y0 [float]: the fiber center y in pixels
        x0 [float]: the fiber center x in pixels
        radius [float]: the fiber radius in pixels
    """
    if method is None:
        if image_obj.camera == 'ff':
            method = 'gaussian'
        elif image_obj.camera == 'nf':
            method = 'edge'
        elif image_obj.camera == 'in':
            method = 'edge'
        else:
            raise RuntimeError('Method or camera type must be declared to get image data')
    y0, x0 = image_obj.getFiberCenter(method=method, tol=1,
                                      test_range=10, show_image=False)
    radius = image_obj.getFiberRadius()        
    image_array = image_obj.getImageArray()
    return image_array, y0, x0, radius

#=============================================================================#
#==== Scrambling Gain Functions ==============================================#
#=============================================================================#

def scramblingGain(in_obj, out_obj, method='edge'):
    """Calculates the scrambling gain for fiber input and output images

    Args:
        in_obj [list(ImageAnalysis)]: list of the ImageAnalysis input objects
        out_obj [list(ImageAnalysis)]: list of the ImageAnalysis output objects

    Returns:
        scrambling_gain [float]: approximate scrambling gain
        input_x [list]
        input_y [list]
        output_x [list]
        output_y [list]
    """
    in_centroid_1 = in_obj_1.getFiberCentroid()
    in_diameter_1 = in_obj_1.getFiberDiameter(method=method)

    out_centroid_1 = out_obj_1.getFiberCentroid()
    out_diameter_1 = out_obj_1.getFiberDiameter(method=method)

    if in_obj_2 is None or out_obj_2 is None:
        in_centroid_2 = in_obj_1.getFiberCenter(method=method)
        out_centroid_2 = out_obj_1.getFiberCenter(method=method)
        in_diameter = in_diameter_1
        out_diameter = out_diameter_1

    else:
        in_centroid_2 = in_obj_2.getFiberCentroid()
        out_centroid_2 = out_obj_2.getFiberCentroid()
        in_diameter = (in_diameter_1 + in_obj_2.getFiberDiameter(method=method)) / 2.0
        out_diameter = (out_diameter_1 + in_obj_2.getFiberDiameter(method=method)) / 2.0


    delta_D_in = np.sqrt((in_centroid_1[1] - in_centroid_2[1])**2 + (in_centroid_1[0] - in_centroid_2[0])**2)
    delta_D_out = np.sqrt((out_centroid_1[1] - out_centroid_2[1])**2 + (out_centroid_1[0] - out_centroid_2[0])**2)

    scramblingGain = (delta_D_in / in_diameter) / (delta_D_out / out_diameter)

    return scramblingGain

#=============================================================================#
#==== Focal Ratio Degradation Functions ======================================#
#=============================================================================#

def FRD(ff_obj, input_f_ratio=-1.0, f_lim=(2.4, 10.0), res=0.1):
    """Calculates the encircled energy for various f ratios

    Args:
        ff_obj [ImageAnalysis]: the image object on which the FRD is calculated
        input_f_ratio [float]: the fiber input f ratio
        f_lim [(float, float)]: the limits of the f ratio calculations
        res [float]: the spacing between each f ratio when calculating
            encircld energy

    Returns:
        f_ratios [list]: list of the f ratios used to calculate encircled
            energy
        encircled_energy [list]: list of the encircled energy at each given
            f ratio
        energy_loss [float]: the loss of energy when the output f ratio equals
            the input f ratio given as a percent
        output_f_ratio [float]: the approximate f ratio inside which 95% of 
            the total encircled energy is included
    """
    center_y, center_x = ff_obj.getFiberCentroid()

    f_ratios = np.arange(f_lim[0], f_lim[1], res)
    energy_loss = None
    output_f_ratio = None
    encircled_energy = []
    for f in f_ratios:
        radius = _f_ratio_to_radius(f, ff_obj.magnification)
        isolated_circle = ff_obj.isolateCircle(ff_obj.getImage(),
                                               center_x,
                                               center_y,
                                               radius)
        iso_circ_sum = isolated_circle.sum()
        encircled_energy.append(iso_circ_sum)
        if abs(f - input_f_ratio) < res / 2.0:
            energy_loss = 100 * (1 - iso_circ_sum / encircled_energy[0])
        if iso_circ_sum / encircled_energy[0] >= 0.95:
            output_f_ratio = f

    encircled_energy = list(np.array(encircled_energy) / encircled_energy[0])

    return f_ratios, encircled_energy, energy_loss, output_f_ratio

def _f_ratio_to_radius(f_ratio, magnification):
    """Converts an f ratio to a radius in units of pixels
    
    Args:
        f_ratio [float]: the f ratio to be converted
        magnification [float]: the magnification of the camera

    Returns:
        radius [float]: in units of pixels
    """
    return magnification * 4.0 * 25400 / f_ratio / 2.0 / 3.45

#=============================================================================#
#==== Modal Noise Functions ===================================================#
#=============================================================================#

def modalNoise(self, image_obj=None, method='fft', output='array', radius_factor=0.95, deg=4):
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
        return self._modalNoiseTophat(image_obj, output, radius_factor)
    elif method in 'fft' or method in 'fourier':
        return self._modalNoiseFFT(image_obj, output, radius_factor)
    elif method in 'polynomial':
        return self._modalNoisePolynomial(image_obj, output, radius_factor, deg=deg)
    elif method in 'gaussian':
        return self._modalNoiseGaussian(image_obj, output, radius_factor)
    elif method in 'gradient':
        return self._modalNoiseGradient(image_obj, output, radius_factor)
    elif method in 'contrast':
        return self._modalNoiseContrast(image_obj, radius_factor)
    elif method in 'gini':
        return self._modalNoiseGini(image_obj, radius_factor)
    elif method in 'entropy':
        return self._modalNoiseEntropy(image_obj, radius_factor)
    else:
        raise ValueError('Incorrect string for method type')

def _modalNoiseFFT(self, image_obj=None, output='array', radius_factor=1.05):
    """Finds modal noise of image using the image's power spectrum
    
    Args:
        image_obj (optional): ImageAnalysis object to be used
        output: chooses whether to output as 'array' or 'parameter'
        radius_factor (optional): fraction of the radius outside which
            the array is padded with zeros

    Returns:
        array: (fft_array, freq_array) where fft_array is the normalized
            azimuthally averaged power spectrum and freq_array are the
            respective frequencies in 1/um
        parameter: finds the Gini coefficient for the 2D power spectrum
    """
    if image_obj is None:
        image_obj = self.image_obj

    image_array, y0, x0, radius = self.getImageData(image_obj)
    height, width = image_array.shape

    if image_obj.camera == 'nf':
        image_array, x0, y0 = self.cropImage(image_array, x0, y0, radius*radius_factor)
        image_array = self.isolateCircle(image_array, x0, y0, radius*radius_factor)

    elif image_obj.camera == 'ff':
        image_array, x0, y0 = self.cropImage(image_array, x0, y0, min(x0, y0, width-x0, height-y0))

    image_array = self.applyWindow(image_array)
    height, width = image_array.shape        

    fft_length = 2500 #8 * min(height, width)
    fft_array = np.fft.fftshift(np.abs(np.fft.fft2(image_array, s=(fft_length, fft_length), norm='ortho')))
    fx0 = fft_length/2
    fy0 = fft_length/2

    max_freq = fft_length/2     

    if len(output) < 3:
        raise ValueError('Incorrect output string')

    elif output in 'array':
        bin_width = 1
        list_len = int(max_freq / bin_width) + 1

        fft_list = np.zeros(list_len).astype('float64')
        weight_list = np.zeros(list_len).astype('float64')
        freq_list = bin_width * np.arange(list_len).astype('float64')

        bottom_right = fft_array[fy0:fy0+max_freq, fx0:fx0+max_freq]
        bottom_left = fft_array[fy0:fy0+max_freq, fx0-max_freq+1:fx0+1][:, ::-1]
        top_left = fft_array[fy0-max_freq+1:fy0+1, fx0-max_freq+1:fx0+1][::-1, ::-1]
        top_right = fft_array[fy0-max_freq+1:fy0+1, fx0:fx0+max_freq][::-1, :]

        fft_array = (bottom_right + bottom_left + top_left + top_right) / 4.0
        #self.showImageArray(np.log(fft_array))

        for i in xrange(max_freq):
            for j in xrange(i+1):
                freq = np.sqrt(i**2 + j**2)
                if freq <= max_freq:
                    fft_list[int(freq/bin_width)] += fft_array[j, i] + fft_array[i, j]
                    weight_list[int(freq/bin_width)] += 2.0

        # Remove bins with nothing in them
        mask = (weight_list > 0.0).astype('bool')
        weight_list = weight_list[mask]
        freq_list = freq_list[mask]
        fft_list = fft_list[mask] / weight_list # Average out

        fft_list /= fft_list.sum() # Normalize
        freq_list /= fft_length * image_obj.pixel_size / image_obj.magnification # Get frequencies in 1/um

        self.plotFFT([freq_list], [fft_list], labels=['Modal Noise Power Spectrum'])

        return fft_list, freq_list

    elif output in 'parameter':
        return self._gini_coefficient(intensityArray(fft_array, fx0, fy0, max_freq))

    else:
        ValueError('Incorrect output string')

def _modalNoiseTophat(self, image_obj=None, output='array', radius_factor=0.95):
    """Finds modal noise of image assumed to be a tophat
    
    Modal noise is defined as the variance across the fiber face normalized
    by the mean across the fiber face

    Args:
        image_obj (optional): ImageAnalysis object to be used
        output: chooses whether to output as 'array' or 'parameter'
        radius_factor (optional): fraction of the radius inside which the
            modal noise calculation will be made

    Returns:
        array: (circle_fit, image_array) 2D numpy arrays of the best fit
            tophat and the fiber image zoomed in to relevant area
        parameter: STDEV / MEAN for the intensities inside the fiber face
    """
    if image_obj is None:
        image_obj = self.image_obj

    image_array, y0, x0, radius = self.getImageData(image_obj)
   
    intensity_array = intensityArray(image_array, x0, y0, radius*radius_factor)
    
    if len(output) < 3:
        raise ValueError('Incorrect output string')

    elif output in 'array': 
        circle_fit = intensity_array.mean() * self.circleArray(self.getMeshGrid(image_array),
                                                               x0, y0, radius, res=10)
        image_array, x0_new, y0_new = self.cropImage(image_array, x0, y0, radius*self.plot_buffer)
        circle_fit = self.cropImage(circle_fit, x0, y0, radius*self.plot_buffer)[0]
        self.showImageArray(image_array)
        self.showImageArray(circle_fit)

        self.plotOverlaidCrossSections(image_array, circle_fit, y0_new, x0_new)
        return circle_fit, image_array

    elif output in 'parameter':
        return intensity_array.std() / intensity_array.mean()

    else:
        raise ValueError('Incorrect output string')

def _modalNoiseGradient(self, image_obj=None, output='array', radius_factor=0.95):
    """Finds modal noise of image using the image gradient

    Args:
        image_obj (optional): ImageAnalysis object to be used
        output: chooses whether to output as 'array' or 'parameter'
        radius_factor (optional): fraction of the radius inside which the
            modal noise calculation will be made

    Returns:
        array: (gradient_array, image_array) 2D numpy array representing magnitude of the gradient
        parameter: STDEV / MEAN for the gradient in the fiber image
    """
    if image_obj is None:
        image_obj = self.image_obj

    image_array, y0, x0, radius = self.getImageData(image_obj)
    image_array, x0, y0 = self.cropImage(image_array, x0, y0, radius*self.plot_buffer)

    gradient_y, gradient_x = np.gradient(image_array)
    gradient_array = np.sqrt(gradient_x**2 + gradient_y**2)

    if len(output) < 3:
        raise ValueError('Incorrect output string')

    elif output in 'array':
        self.showImageArray(gradient_array)
        self.plotOverlaidCrossSections(image_array, gradient_array, y0, x0)
        return gradient_array, image_array

    elif output in 'parameter':
        intensity_array = intensityArray(gradient_array, x0, y0, radius*radius_factor)
        image_intensity_array = intensityArray(image_array, x0, y0, radius*radius_factor)
        return intensity_array.std() / image_intensity_array.mean()

    else:
        ValueError('Incorrect output string')

def _modalNoisePolynomial(self, image_obj=None, output='array', radius_factor=0.95, deg=4):
    """Finds modal noise of image using polynomial fit
    
    Crops image exactly around the circumference of the circle and fits a
    polynomial to the 2D image. Modal noise is then defined as the STDEV
    in the difference between the original image and the polynomial fit
    normalized by the mean value inside the fiber face

    Args:
        image_obj (optional): ImageAnalysis object to be used            
        output: chooses whether to output as 'array' or 'parameter'
        radius_factor (optional): fraction of the radius inside which the
            modal noise calculation will be made

    Returns:
        array: (poly_fit, image_array) 2D numpy arrays of the best fit
            polynomial and the fiber image zoomed in to relevant area
        parameter: STDEV for the difference between the fiber image and
            poly_fit divided by the mean of the fiber image intensities
    """
    if image_obj is None:
        image_obj = self.image_obj

    image_array, y0, x0, radius = self.getImageData(image_obj)
    radius *= radius_factor / np.sqrt(2)

    image_array, x0, y0 = self.cropImage(image_array, x0, y0, radius)

    poly_fit = self.getPolynomialFit(image_array, deg=deg, x0=x0, y0=y0)

    if len(output) < 3:
        raise ValueError('Incorrect output string')

    elif output in 'array':
        self.plotOverlaidCrossSections(image_array, poly_fit, radius, radius)
        return poly_fit, image_array

    elif output in 'parameter':
        diff_array = image_array - poly_fit

        intensity_array = intensityArray(diff_array, x0, y0, radius * np.sqrt(2))
        image_intensity_array = intensityArray(image_array, x0, y0, radius * np.sqrt(2))

        return intensity_array.std() / image_intensity_array.mean()

    else:
        raise ValueError('Incorrect output string')

def _modalNoiseGaussian(self, image_obj=None, output='array', radius_factor=0.95):
    """Finds modal noise of image using gaussian fit
    
    Crops image exactly around the circumference of the circle and fits a
    gaussian to the 2D image. Modal noise is then defined as the STDEV
    in the difference between the original image and the gaussian fit
    normalized by the mean value inside the fiber face

    Args:
        image_obj (optional): ImageAnalysis object to be used            
        output: chooses whether to output as 'array' or 'parameter'
        radius_factor (optional): fraction of the radius inside which the
            modal noise calculation will be made

    Returns:
        array: (gauss_fit, image_array) 2D numpy arrays of the best fit
            polynomial and the fiber image zoomed in to relevant area
        parameter: STDEV for the difference between the fiber image and
            gauss_fit divided by the mean of the fiber image intensities
    """
    if image_obj is None:
        image_obj = self.image_obj

    image_array, y0, x0, radius = self.getImageData(image_obj)
    radius *= radius_factor / np.sqrt(2)

    image_array, x0, y0 = self.cropImage(image_array, x0, y0, radius)

    gauss_fit = self.getGaussianFit(image_array)

    if len(output) < 3:
        raise ValueError('Incorrect output string')

    elif output in 'array':
        self.showImageArray(gauss_fit)
        self.plotOverlaidCrossSections(image_array, gauss_fit, radius, radius)
        return gauss_fit, image_array

    elif output in 'parameter':
        diff_array = image_array - gauss_fit

        intensity_array = intensityArray(diff_array, x0, y0, radius)
        image_intensity_array = intensityArray(image_array, x0, y0, radius)

        return intensity_array.std() / image_intensity_array.mean()

    else:
        raise ValueError('Incorrect output string')

def _modalNoiseGini(self, image_obj=None, radius_factor=0.95):
    """Find modal noise of image using Gini coefficient

    Args:
        image_obj (optional): ImageAnalysis object to be used
        radius_factor (optional): fraction of the radius inside which the
            modal noise calculation will be made

    Returns:
        array: 2D numpy array zoomed into the fiber face
        parameter: Gini coefficient for intensities inside fiber face
    """
    if image_obj is None:
        image_obj = self.image_obj

    image_array, y0, x0, radius = self.getImageData(image_obj)
    radius *= radius_factor        

    intensity_array = intensityArray(image_array, x0, y0, radius)

    return self._gini_coefficient(intensity_array)

def _gini_coefficient(self, test_array):
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

def _modalNoiseContrast(self, image_obj=None, radius_factor=0.95):
    """Finds modal noise of image using Michelson contrast
    
    Modal noise is defined as (I_max - I_min) / (I_max + I_min)

    Args:
        image_obj (optional): ImageAnalysis object to be used
        radius_factor (optional): fraction of the radius inside which the
            modal noise calculation will be made

    Returns:
        parameter: (I_max - I_min) / (I_max + I_min) for intensities inside
            fiber face
    """
    if image_obj is None:
        image_obj = self.image_obj

    image_array, y0, x0, radius = self.getImageData(image_obj)
    radius *= radius_factor        
       
    intensity_array = intensityArray(image_array, x0, y0, radius)

    return (intensity_array.max() - intensity_array.min()) / (intensity_array.max() + intensity_array.min())     

def _modalNoiseEntropy(self, image_obj=None, radius_factor=0.95):
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

    intensity_array = intensityArray(image_array, x0, y0, radius)
    intensity_array = intensity_array / intensity_array.sum()
    return (-intensity_array * np.log10(intensity_array)).sum()

#=============================================================================#
#==== Modal Noise Test =======================================================#
#=============================================================================#

if __name__ == '__main__':
    from ImageAnalysis import ImageAnalysis
    from Calibration import Calibration
    import os as os
    import matplotlib.pyplot as plt
    from copy import deepcopy
    plt.rc('font', size=16, family='sans-serif')
    plt.rc('xtick', labelsize=14)
    plt.rc('ytick', labelsize=14)
    plt.rc('lines', lw=2)

    base_folder = '2016-07-26/'
    ambient_folder = base_folder + 'ambient/600um/'
    dark_folder = base_folder + 'dark/'
    agitated_folder = base_folder + 'single/600um/agitated/'
    unagitated_folder = base_folder + 'single/600um/unagitated/'
    ext = '.fit'

    nf = {}
    ff = {}

    nf['calibration'] = Calibration([dark_folder + 'nf_dark_' + str(i).zfill(3) + ext for i in xrange(10)],
                                    None,
                                    [ambient_folder + 'nf_ambient_' + str(i).zfill(3) + '_0.1' + ext for i in xrange(10)])
    print 'NF calibration initialized'
    ff['calibration'] = Calibration([dark_folder + 'ff_dark_' + str(i).zfill(3) + ext for i in xrange(10)],
                                    None,
                                    [ambient_folder + 'ff_ambient_' + str(i).zfill(3) + '_0.1' + ext for i in xrange(10)])
    print 'FF calibration initialized'

    empty_data = {'images': [], 'fft': [], 'freq': []}
    for test in ['agitated', 'unagitated', 'baseline']:
        nf[test] = deepcopy(empty_data)
        ff[test] = deepcopy(empty_data)

    image_range = xrange(20)
    nf['agitated']['images'] = [agitated_folder + 'nf_agitated_' + str(i).zfill(3) + ext for i in image_range]
    nf['unagitated']['images'] = [unagitated_folder + 'nf_unagitated_' + str(i).zfill(3) + ext for i in image_range]

    ff['agitated']['images'] = [agitated_folder + 'ff_agitated_' + str(i).zfill(3) + ext for i in image_range]
    ff['unagitated']['images'] = [unagitated_folder + 'ff_unagitated_' + str(i).zfill(3) + ext for i in image_range]

    for test in ['agitated', 'unagitated']:
        nf[test]['obj'] = ImageAnalysis(nf[test]['images'], nf['calibration'])

    nf_test_obj = nf['agitated']['obj']
    nf['baseline']['obj'] = ImageAnalysis(nf_test_obj.getTophatFit(),
                                          pixel_size=nf_test_obj.pixel_size,
                                          camera='nf',
                                          threshold = 0.1)

    for test in ['agitated', 'unagitated']:
        ff[test]['obj'] = ImageAnalysis(ff[test]['images'], ff['calibration'])
    ff_test_obj = ff['agitated']['obj']
    ff['baseline']['obj'] = ImageAnalysis(ff_test_obj.getGaussianFit(),
                                          pixel_size=ff_test_obj.pixel_size,
                                          camera='ff')

    print

    for cam in [nf, ff]:
        for test in ['agitated', 'unagitated', 'baseline']:
            print test + ':'
            cam[test]['fft'], cam[test]['freq'] = modalNoise(cam[test]['obj'],
                                                             method='fft',
                                                             output='array',
                                                             radius_factor=1.0)

        NumpyArrayHandler.plotFFT([cam[test]['freq'] for test in ['agitated', 'unagitated', 'baseline']],
                                  [cam[test]['fft'] for test in ['agitated', 'unagitated', 'baseline']],
                                  ['Agitated laser', 'Unagitated laser', 'Baseline'],
                                  cam[test]['obj'].camera.upper() + ' Modal Noise Comparison (600um Fiber)')

