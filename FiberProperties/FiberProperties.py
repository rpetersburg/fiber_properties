"""FiberProperties.py was written by Ryan Petersburg for use with fiber
characterization for the EXtreme PRecision Spectrograph

This module contains functions that calculate the scrambling gain, focal ratio
degradation, and modal noise for image analysis objects
"""
import numpy as np
import matplotlib.pyplot as plt
from NumpyArrayHandler import *

def getImageData(image_obj, method=None):
    """Returns relevant information from an ImageAnalysis object
    
    Args:
        image_obj [ImageAnalysis]: image object to be analyzed

    Returns:
        image_array [ndarray]: the 2D image
        y0 [float]: the fiber center y in pixels
        x0 [float]: the fiber center x in pixels
        radius [float]: the fiber radius in pixels
    """
    if method is None:
        if image_obj.getCamera() == 'ff':
            method = 'gaussian'
        elif image_obj.getCamera() == 'nf':
            method = 'edge'
        elif image_obj.getCamera() == 'in':
            method = 'edge'
        else:
            raise RuntimeError('Method or camera type must be declared to get image data')
    y0, x0 = image_obj.getFiberCenter(method=method, tol=1,
                                      test_range=10, show_image=False)
    radius = image_obj.getFiberRadius(method=method)        
    image_array = image_obj.getImage()
    return image_array, y0, x0, radius

#=============================================================================#
#==== Scrambling Gain Functions ==============================================#
#=============================================================================#

def scramblingGain(in_objs, out_objs, in_method='gaussian', out_method='edge'):
    """Calculates the scrambling gain for fiber input and output images

    Args:
        in_objs [list(ImageAnalysis)]: list of the ImageAnalysis input objects
        out_objs [list(ImageAnalysis)]: list of the ImageAnalysis output objects

    Returns:
        input_x [list]
        input_y [list]
        output_x [list]
        output_y [list]        
        scrambling_gain [float]: approximate scrambling gain
    """
    if not isinstance(in_objs, Iterable) or not isinstance(out_objs, Iterable):
        in_objs = [in_objs]
        out_objs = [out_objs]

    if len(in_objs) != len(out_objs):
        raise RuntimeError('Lists of input and output objects not the same length')

    input_x = []
    input_y = []
    for in_obj in in_objs:
        in_centroid = in_obj.getFiberCentroid(radius_factor=1.05, method=in_method, units='microns')
        in_center = in_obj.getFiberCenter(method=out_method, units='microns')
        in_diameter = in_obj.getFiberDiameter(method=out_method, units='microns')
        input_x.append((in_centroid[1] - in_center[1]) / in_diameter)
        input_y.append((in_centroid[0] - in_center[0]) / in_diameter)

    output_x = []
    output_y = []
    output_diameter = out_objs[0].getFiberDiameter(method=out_method, units='microns')
    for out_obj in out_objs:
        out_centroid = out_obj.getFiberCentroid(radius_factor=1.0, method=out_method, units='microns')
        out_center = out_obj.getFiberCenter(method=out_method, units='microns')
        out_diameter = out_obj.getFiberDiameter(method=out_method, units='microns')
        output_x.append((out_centroid[1] - out_center[1]) / out_diameter)
        output_y.append((out_centroid[0] - out_center[0]) / out_diameter)

    scrambling_gain = None

    return input_x, input_y, output_x, output_y, scrambling_gain

def twoPointScramblingGain(in_centroid_1, in_centroid_2, out_centroid_1, out_centroid_2):
    return

#=============================================================================#
#==== Focal Ratio Degradation Functions ======================================#
#=============================================================================#

def FRD(ff_obj, input_focal_ratio=-1.0, focal_lim=(2.4, 10.0), res=0.1):
    """Calculates the encircled energy for various f ratios

    Args:
        ff_obj [ImageAnalysis]: the image object on which the FRD is calculated
        input_focal_ratio [float]: the fiber input f ratio
        focal_lim [(float, float)]: the limits of the f ratio calculations
        res [float]: the spacing between each f ratio when calculating
            encircld energy

    Returns:
        focal_ratios [list]: list of the f ratios used to calculate encircled
            energy
        encircled_energy [list]: list of the encircled energy at each given
            f ratio
        energy_loss [float]: the loss of energy when the output f ratio equals
            the input f ratio given as a percent
        output_focal_ratio [float]: the approximate f ratio inside which 95% of 
            the total encircled energy is included
    """
    center_y, center_x = ff_obj.getFiberCentroid()

    focal_ratios = list(np.arange(focal_lim[0], focal_lim[1] + res, res))
    energy_loss = None
    output_focal_ratio = None
    encircled_energy = []
    for f in focal_ratios:
        radius = _focal_ratio_to_radius(f, ff_obj)
        isolated_circle = isolateCircle(ff_obj.getImage(),
                                        center_x,
                                        center_y,
                                        radius)
        iso_circ_sum = isolated_circle.sum()
        encircled_energy.append(iso_circ_sum)
        if abs(f - input_focal_ratio) < res / 2.0:
            energy_loss = 100 * (1 - iso_circ_sum / encircled_energy[0])
        if iso_circ_sum / encircled_energy[0] >= 0.95:
            output_focal_ratio = f

    encircled_energy = list(np.array(encircled_energy) / encircled_energy[0])

    return focal_ratios, encircled_energy, energy_loss, output_focal_ratio

def _focal_ratio_to_radius(focal_ratio, im_obj):
    """Converts an f ratio to an image radius in units of pixels
    
    Args:
        focal_ratio [float]: the f ratio to be converted
        magnification [float]: the magnification of the camera

    Returns:
        radius [float]: in units of pixels
    """
    return 25400 * (4.0 / focal_ratio) * (im_obj.magnification / im_obj.pixel_size) / 2.0

#=============================================================================#
#==== Modal Noise Functions ===================================================#
#=============================================================================#

def modalNoise(image_obj, method='fft', **kwargs):
    """Finds modal noise of image using specified method and output

    Args:
        image_obj [ImageAnalysis]: the image object being analyzed
        method [{'tophat', 'fft', 'polynomial', 'gaussian', 'gradient',
                 'contrast', 'gini', 'entropy'}]:
            string designating the modal noise method to use
        **kwargs :
            The keyworded arguments to pass to the modal noise method

    Returns:
        modal_noise_parameter: if output is 'parameter'
        modal_noise_arrays: if output is 'arrays'
    """
    if image_obj is None:
        image_obj = image_obj

    if len(method) < 3:
        raise ValueError('Incorrect string for method type')
    elif method in 'tophat':
        return _modalNoiseTophat(image_obj, **kwargs)
    elif method in 'fft' or method in 'fourier':
        return _modalNoiseFFT(image_obj, **kwargs)
    elif method in 'polynomial':
        return _modalNoisePolynomial(image_obj, **kwargs)
    elif method in 'gaussian':
        return _modalNoiseGaussian(image_obj, **kwargs)
    elif method in 'gradient':
        return _modalNoiseGradient(image_obj, **kwargs)
    elif method in 'contrast':
        return _modalNoiseContrast(image_obj, **kwargs)
    elif method in 'gini':
        return _modalNoiseGini(image_obj, **kwargs)
    elif method in 'entropy':
        return _modalNoiseEntropy(image_obj, **kwargs)
    else:
        raise ValueError('Incorrect string for method type')

def _modalNoiseFFT(image_obj, output='array', radius_factor=1.05):
    """Finds modal noise of image using the image's power spectrum
    
    Args:
        image_obj [ImageAnalysis]: image object to analyze
        output [{'array', 'parameter'}, optional]: see below for further info
        radius_factor [number, optional]: fraction of the radius outside which
            the array is padded with zeros

    Returns:
        output='array':
            fft_array [np.ndarray]: normalized, azimuthally averaged power
                spectrum
            freq_array [np.ndarray]: respective frequencies in 1/um
        output='parameter':
            parameter: the Gini coefficient for the 2D power spectrum
    """
    image_array, y0, x0, radius = getImageData(image_obj)
    height, width = image_array.shape

    if image_obj.getCamera() == 'nf':
        image_array, x0, y0 = cropImage(image_array, x0, y0, radius*radius_factor)
        image_array = isolateCircle(image_array, x0, y0, radius*radius_factor)

    elif image_obj.getCamera() == 'ff':
        image_array, x0, y0 = cropImage(image_array, x0, y0, min(x0, y0, width-x0, height-y0))

    image_array = applyWindow(image_array)
    height, width = image_array.shape        

    fft_length = 2500 #8 * min(height, width)
    fft_array = np.fft.fftshift(np.abs(np.fft.fft2(image_array, s=(fft_length, fft_length), norm='ortho')))
    showImageArray(image_obj.getImage())
    showImageArray(np.log(fft_array))
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

        # Take the four quadrants of the FFT and sum them together
        bottom_right = fft_array[fy0:fy0+max_freq, fx0:fx0+max_freq]
        bottom_left = fft_array[fy0:fy0+max_freq, fx0-max_freq+1:fx0+1][:, ::-1]
        top_left = fft_array[fy0-max_freq+1:fy0+1, fx0-max_freq+1:fx0+1][::-1, ::-1]
        top_right = fft_array[fy0-max_freq+1:fy0+1, fx0:fx0+max_freq][::-1, :]

        fft_array = (bottom_right + bottom_left + top_left + top_right) / 4.0

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
        freq_list /= fft_length * image_obj.getPixelSize() / image_obj.getMagnification() # Get frequencies in 1/um

        plotFFT([freq_list], [fft_list], labels=['Modal Noise Power Spectrum'])

        return fft_list, freq_list

    elif output in 'parameter':
        return _gini_coefficient(intensityArray(fft_array, fx0, fy0, max_freq))

    else:
        raise ValueError('Incorrect output string')

def _modalNoiseTophat(image_obj, output='array', radius_factor=0.95):
    """Finds modal noise of image assumed to be a tophat
    
    Modal noise is defined as the variance across the fiber face normalized
    by the mean across the fiber face

    Args:
        image_obj [ImageAnalysis]: image object to analyze
        output [{'array', 'parameter'}, optional]: see below for further info
        radius_factor [number, optional]: fraction of the radius inside which
            the modal noise calculation will be made

    Returns:
        output='array':
            circle_fit [ndarray]: best fit 2D tophat image
            image_array [ndarray]: fiber image zoomed to relevant area
        output='parameter':
            parameter [float]: STDEV / MEAN for the intensities inside the
                fiber face
    """
    image_array, y0, x0, radius = getImageData(image_obj)
   
    intensity_array = intensityArray(image_array, x0, y0, radius*radius_factor)
    
    if len(output) < 3:
        raise ValueError('Incorrect output string')

    elif output in 'array': 
        circle_fit = intensity_array.mean() * circleArray(getMeshGrid(image_array),
                                                               x0, y0, radius, res=10)
        image_array, x0_new, y0_new = cropImage(image_array, x0, y0, radius*plot_buffer)
        circle_fit = cropImage(circle_fit, x0, y0, radius*plot_buffer)[0]
        showImageArray(image_array)
        showImageArray(circle_fit)

        plotOverlaidCrossSections(image_array, circle_fit, y0_new, x0_new)
        return circle_fit, image_array

    elif output in 'parameter':
        return intensity_array.std() / intensity_array.mean()

    else:
        raise ValueError('Incorrect output string')

def _modalNoiseGradient(image_obj, output='array', radius_factor=0.95):
    """Finds modal noise of image using the image gradient

    Args:
        image_obj [ImageAnalysis]: image object to analyze
        output [{'array', 'parameter'}, optional]: see below for further info
        radius_factor [number, optional]: fraction of the radius inside which
            the modal noise calculation will be made

    Returns:
        output='array':
            gradient_array [ndarray]: 2D gradient magnitude image
            image_array [ndarray]: fiber image zoomed to relevant area
        output='parameter'
            parameter [float]: STDEV / MEAN for the gradient in the fiber image
    """
    if image_obj is None:
        image_obj = image_obj

    image_array, y0, x0, radius = getImageData(image_obj)
    image_array, x0, y0 = cropImage(image_array, x0, y0, radius*plot_buffer)

    gradient_y, gradient_x = np.gradient(image_array)
    gradient_array = np.sqrt(gradient_x**2 + gradient_y**2)

    if len(output) < 3:
        raise ValueError('Incorrect output string')

    elif output in 'array':
        showImageArray(gradient_array)
        plotOverlaidCrossSections(image_array, gradient_array, y0, x0)
        return gradient_array, image_array

    elif output in 'parameter':
        intensity_array = intensityArray(gradient_array, x0, y0, radius*radius_factor)
        image_intensity_array = intensityArray(image_array, x0, y0, radius*radius_factor)
        return intensity_array.std() / image_intensity_array.mean()

    else:
        ValueError('Incorrect output string')

def _modalNoisePolynomial(image_obj, output='array', radius_factor=0.95, deg=4):
    """Finds modal noise of image using polynomial fit
    
    Crops image exactly around the circumference of the circle and fits a
    polynomial to the 2D image. Modal noise is then defined as the STDEV
    in the difference between the original image and the polynomial fit
    normalized by the mean value inside the fiber face

    Args:
        image_obj [ImageAnalysis]: image object to analyze
        output [{'array', 'parameter'}, optional]: see below for further info
        radius_factor [number, optional]: fraction of the radius inside which
            the modal noise calculation will be made

    Returns:
        output='array':
            poly_fit [ndarray]: best fit 2D polynomial image
            image_array [ndarray]: fiber image zoomed to relevant area
        output='parameter':
            parameter [float]: STDEV for the difference between the fiber image
                and poly_fit divided by the mean of the fiber image intensities
    """
    if image_obj is None:
        image_obj = image_obj

    image_array, y0, x0, radius = getImageData(image_obj)
    radius *= radius_factor / np.sqrt(2)

    image_array, x0, y0 = cropImage(image_array, x0, y0, radius)

    poly_fit = getPolynomialFit(image_array, deg=deg, x0=x0, y0=y0)

    if len(output) < 3:
        raise ValueError('Incorrect output string')

    elif output in 'array':
        plotOverlaidCrossSections(image_array, poly_fit, radius, radius)
        return poly_fit, image_array

    elif output in 'parameter':
        diff_array = image_array - poly_fit

        intensity_array = intensityArray(diff_array, x0, y0, radius * np.sqrt(2))
        image_intensity_array = intensityArray(image_array, x0, y0, radius * np.sqrt(2))

        return intensity_array.std() / image_intensity_array.mean()

    else:
        raise ValueError('Incorrect output string')

def _modalNoiseGaussian(image_obj, output='array', radius_factor=0.95):
    """Finds modal noise of image using a gaussian fit
    
    Crops image exactly around the circumference of the circle and fits a
    gaussian to the 2D image. Modal noise is then defined as the STDEV
    in the difference between the original image and the gaussian fit
    normalized by the mean value inside the fiber face

    Args:
        image_obj [ImageAnalysis]: image object to analyze
        output [{'array', 'parameter'}, optional]: see below for further info
        radius_factor [number, optional]: fraction of the radius inside which
            the modal noise calculation will be made

    Returns:
        output='array':
            gauss_fit [ndarray]: best fit 2D gaussian image
            image_array [ndarray]: fiber image zoomed to relevant area
        output='parameter'
            parameter [float]: STDEV for the difference between the fiber image
                and gauss_fit divided by the mean of the fiber image
                intensities
    """
    if image_obj is None:
        image_obj = image_obj

    image_array, y0, x0, radius = getImageData(image_obj)
    radius *= radius_factor / np.sqrt(2)

    image_array, x0, y0 = cropImage(image_array, x0, y0, radius)

    gauss_fit = getGaussianFit(image_array)

    if len(output) < 3:
        raise ValueError('Incorrect output string')

    elif output in 'array':
        showImageArray(gauss_fit)
        plotOverlaidCrossSections(image_array, gauss_fit, radius, radius)
        return gauss_fit, image_array

    elif output in 'parameter':
        diff_array = image_array - gauss_fit

        intensity_array = intensityArray(diff_array, x0, y0, radius)
        image_intensity_array = intensityArray(image_array, x0, y0, radius)

        return intensity_array.std() / image_intensity_array.mean()

    else:
        raise ValueError('Incorrect output string')

def _modalNoiseGini(image_obj, radius_factor=0.95):
    """Find modal noise of image using Gini coefficient

    Args:
        image_obj [ImageAnalysis]: image object to analyze
        radius_factor [number, optional]: fraction of the radius inside which
            the modal noise calculation will be made

    Returns:
        parameter [float]: Gini coefficient for intensities inside fiber face
    """
    if image_obj is None:
        image_obj = image_obj

    image_array, y0, x0, radius = getImageData(image_obj)
    radius *= radius_factor        

    intensity_array = intensityArray(image_array, x0, y0, radius)

    return _gini_coefficient(intensity_array)

def _gini_coefficient(test_array):
    """Finds gini coefficient for intensities in given array

    Args:
        test_array [ndarray]: arbitrary-dimension numpy array of values

    Returns:
        gini_coefficient [float]
    """
    test_array = test_array.ravel()
    gini_coeff = 0.0
    for i in xrange(test_array.size):
        gini_coeff += np.abs(test_array[i] - test_array).sum()
    norm = 2.0 * test_array.sum() * test_array.size

    return gini_coeff / norm

def _modalNoiseContrast(image_obj, radius_factor=0.95):
    """Finds modal noise of image using Michelson contrast
    
    Modal noise is defined as (I_max - I_min) / (I_max + I_min)

    Args:
        image_obj [ImageAnalysis]: image object to analyze
        radius_factor [number, optional]: fraction of the radius inside which
            the modal noise calculation will be made

    Returns:
        parameter [float]: (I_max - I_min) / (I_max + I_min) for intensities
            inside fiber face
    """
    if image_obj is None:
        image_obj = image_obj

    image_array, y0, x0, radius = getImageData(image_obj)
    radius *= radius_factor        
       
    intensity_array = intensityArray(image_array, x0, y0, radius)

    return (intensity_array.max() - intensity_array.min()) / (intensity_array.max() + intensity_array.min())     

def _modalNoiseEntropy(image_obj, radius_factor=0.95):
    """Find modal noise of image using hartley entropy

    Args:
        image_obj [ImageAnalysis]: image object to analyze
        radius_factor [number, optional]: fraction of the radius inside which
            the modal noise calculation will be made

    Returns:
        parameter [float]: hartley entropy for intensities inside fiber face
    """
    if image_obj is None:
        image_obj = image_obj

    image_array, y0, x0, radius = getImageData(image_obj)
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
    plt.rc('font', size=14, family='sans-serif')
    plt.rc('xtick', labelsize=14)
    plt.rc('ytick', labelsize=14)
    plt.rc('lines', lw=2)

    base_folder = '../data/Modal Noise Measurements/2016-07-26/'
    ambient_folder = base_folder + 'ambient/600um/'
    dark_folder = base_folder + 'dark/'
    agitated_folder = base_folder + 'images/600um/agitated/'
    unagitated_folder = base_folder + 'images/600um/unagitated/'
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
                                          pixel_size=nf_test_obj.getPixelSize(),
                                          threshold = 0.1,
                                          camera='nf')

    for test in ['agitated', 'unagitated']:
        ff[test]['obj'] = ImageAnalysis(ff[test]['images'], ff['calibration'])
    ff_test_obj = ff['agitated']['obj']
    ff['baseline']['obj'] = ImageAnalysis(ff_test_obj.getGaussianFit(),
                                          pixel_size=ff_test_obj.getPixelSize(),
                                          magnification=1,
                                          camera='ff')

    print

    for cam in [nf, ff]:
        for test in ['agitated', 'unagitated', 'baseline']:
            print test + ':'
            cam[test]['fft'], cam[test]['freq'] = modalNoise(cam[test]['obj'], method='fft', output='array', radius_factor=1.0)

        plotFFT([cam[test]['freq'] for test in ['agitated', 'unagitated', 'baseline']],
                [cam[test]['fft'] for test in ['agitated', 'unagitated', 'baseline']],
                labels=['Agitated laser', 'Unagitated laser', 'Baseline'],
                title=cam[test]['obj'].getCamera().upper() + ' Modal Noise Comparison (600um Fiber)')

