"""modal_noise.py was written by Ryan Petersburg for use with fiber
characterization for the EXtreme PRecision Spectrograph

This module contains functions that calculate the modal noise for
images taken with the FCS contained in FiberImage objects
"""
from __future__ import division
import numpy as np
from .numpy_array_handler import (crop_image, isolate_circle, apply_window,
                                  mesh_grid_from_array, intensity_array,
                                  filter_image)
from .plotting import (plot_image, plot_fft, show_plots, plot_cross_sections,
                       show_image, plot_overlaid_cross_sections, plot_dot)
from .containers import FFTInfo, Pixel

def modal_noise(image_obj, method='fft', **kwargs):
    """Finds modal noise of image using specified method and output

    Args
    ----
    image_obj : FiberImage
        the image object being analyzed
    method : {'tophat', 'fft', 'polynomial', 'gaussian', 'gradient',
             'contrast', 'gini', 'entropy'}
        string designating the modal noise method to use
    **kwargs : dict
        The keyworded arguments to pass to the modal noise method

    Returns
    -------
    See each individual method (most return float)
    """
    if len(method) < 3:
        raise ValueError('String too short for modal noise method')
    elif method in 'fft' or method in 'fourier':
        return _modal_noise_fft(image_obj, **kwargs)
    elif method in 'filtered':
        return _modal_noise_filter(image_obj, **kwargs)
    elif method in 'tophat':
        return _modal_noise_tophat(image_obj, **kwargs)
    elif method in 'contrast':
        return _modal_noise_contrast(image_obj, **kwargs)
    elif method in 'polynomial':
        return _modal_noise_polynomial(image_obj, **kwargs)
    elif method in 'gaussian':
        return _modal_noise_gaussian(image_obj, **kwargs)
    elif method in 'gradient':
        return _modal_noise_gradient(image_obj, **kwargs)
    elif method in 'gini':
        return _modal_noise_gini(image_obj, **kwargs)
    elif method in 'entropy':
        return _modal_noise_entropy(image_obj, **kwargs)
    else:
        raise ValueError('Incorrect string for modal noise method')

def baseline_image(image_obj, kernel_size=None, stdev=0.01, num_images=10, **kwargs):
    """Return a numpy array of a baseline modal noise image

    Args
    ----
    im_obj : FiberImage
        image object to use for baseline

    Returns
    -------
    baseline_image : 2D np.ndarray
        heavily filtered image with added normal distribution noise
    """
    if kernel_size is None:
        kernel_size = 101
    image, center, radius = _get_image_data(image_obj, **kwargs)

    zero_fill = False
    if kernel_size > int((min(*image.shape) - 2*radius)):
        zero_fill = True # Prevents edge effects due to large filters

    image_crop = crop_image(image, center, radius + kernel_size, False)

    perfect_image = filter_image(image_crop, kernel_size=kernel_size, zero_fill=zero_fill)
    perfect_image *= (perfect_image > 0.0).astype('float64')

    baseline_image = np.zeros_like(perfect_image)
    for i in xrange(num_images):
        # baseline_image += np.random.poisson(perfect_image) / 10.0
        baseline_image += np.random.normal(perfect_image, np.sqrt(perfect_image) + stdev) / num_images
    return baseline_image

def _get_image_data(image_obj, fiber_method=None, **kwargs):
    """Returns relevant information from a FiberImage object

    Args
    ----
    image_obj : FiberImage
        Image object to be analyzed
    fiber_method : str, optional
        method to use when calculating center and radius of fiber face

    Returns
    -------
    image : ndarray
        the 2D image
    y0 : float
        the fiber center y
    x0 : float
        the fiber center x
    radius : float
        the fiber radius
    """
    center, diameter = image_obj.get_fiber_data(method=fiber_method,
                                                units='pixels',
                                                **kwargs)
    radius = diameter / 2.0
    image = image_obj.get_image()
    return image, center, radius

def _get_radius_factor(radius, camera=None):
    """Returns the radius factor for 30 pixels inside circumference"""
    if camera == 'ff':
        return 1.0
    return 1 - 30 / radius

def _modal_noise_fft(image_obj, radius_factor=None, show_image=False,
                     **kwargs):
    """Finds modal noise of image using the image's power spectrum

    Args
    ----
    image_obj : FiberImage
        image object to analyze
    show_image : bool, optional
        whether or not to show images of the modal noise analysis
    radius_factor : number, optional
        fraction of the radius outside which the array is padded with zeros

    Returns
    -------
    fft_info : FFTInfo
        se containers.py fr attribute information
    """
    image, center, radius = _get_image_data(image_obj, **kwargs)
    if radius_factor is None:
        radius_factor = _get_radius_factor(radius, image_obj.camera)
    height, width = image.shape

    if image_obj.get_camera() == 'nf':
        image, center = crop_image(image, center, radius*radius_factor)
        image = isolate_circle(image, center, radius*radius_factor)

    elif image_obj.get_camera() == 'ff':
        image, center = crop_image(image, center, min(center.x, center.y,
                                                      width-center.x,
                                                      height-center.y))

    if show_image:
        plot_image(image)

    image = apply_window(image)
    if show_image:
        plot_image(image)
    height, width = image.shape

    fft_length = 2500 # Chosen to get good resolution in decent time
    fft_array = np.fft.fftshift(np.abs(np.fft.fft2(image, s=(fft_length,
                                                             fft_length),
                                                   norm='ortho')))
    if show_image:
        plot_image(np.log(fft_array))

    fx0 = fft_length//2
    fy0 = fft_length//2

    max_freq = fft_length//2

    bin_width = 1
    list_len = max_freq//bin_width + 1

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

    # Normalize
    fft_list /= fft_list.sum()
    # Get Frequencies in 1/um
    freq_list /= image_obj.convert_pixels_to_units(fft_length, 'microns')

    if show_image:
        max_wavelength = image_obj.get_fiber_radius(method='edge', units='microns')
        plot_fft(FFTInfo(np.array(fft_list), np.array(freq_list)),
                 labels=[],
                 max_wavelength=max_wavelength)
        show_plots()

    return FFTInfo(np.array(fft_list), np.array(freq_list))

def _modal_noise_filter(image_obj, kernel_size=None, show_image=False,
                        radius_factor=None, fiber_shape='circle', **kwargs):
    """Finds SNR of image using a median filter comparison

    Find the difference between the image and the median filtered image. Take
    the standard deviation of this difference inside the fiber face as the
    noise divided into the mean across the fiber face

    Args
    ----
    image_obj : FiberImage
        image object to analyze
    kernel_size : odd int, optional
        kernel side length to use for median filter
    show_image : bool, optional
        whether or not to show images of the modal noise analysis
    radius_factor : float, optional
        fraction of the radius inside which the modal noise is calculated

    Returns
    -------
    signal_to_noise : float
        median of intensities inside the fiber face divided by the STDEV of
        difference between these intensities and a median filtered image
    """
    image, center, radius = _get_image_data(image_obj, **kwargs)

    if radius_factor is None:
        radius_factor = _get_radius_factor(radius, image_obj.camera)
    if kernel_size is None:
        kernel_size = 101

    image, new_center = crop_image(image, center, radius)
    filtered_image = filter_image(image, kernel_size, square=False)
    diff_image = image - filtered_image

    if 'rect' in fiber_shape and image_obj.camera != 'ff':
        mask = image_obj.get_threshold_mask()
        mask = crop_image(mask, center, radius, False).astype('bool')
  
        filt_inten_array = filtered_image[mask]
        diff_inten_array = diff_image[mask]

    else:
        filt_inten_array = intensity_array(filtered_image, new_center,
                                           radius*radius_factor)
        diff_inten_array = intensity_array(diff_image, new_center,
                                           radius*radius_factor)
    if show_image:
        plot_image(filtered_image)
        plot_image(diff_image)
        # temp_image = crop_image(image, center, radius*radius_factor)[0]
        # temp_filtered_image, center = crop_image(filtered_image, center,
        #                                          radius*radius_factor)
        # plot_overlaid_cross_sections(temp_image, temp_filtered_image, center)
        plot_overlaid_cross_sections(image, filtered_image, new_center)
        plot_cross_sections(diff_image, new_center)
        show_plots()

    if image_obj.camera == 'ff':
        return np.max(filt_inten_array) / np.std(diff_inten_array)
    return np.median(filt_inten_array) / np.std(diff_inten_array)

def _modal_noise_tophat(image_obj, show_image=False, radius_factor=None, **kwargs):
    """Finds SNR of image assumed to be a tophat

    SNR is defined as the median across the fiber face divided by the standard
    deviation across the fiber face

    Args
    ----
    image_obj : FiberImage
        image object to analyze
    show_image : bool, optional
        whether or not to show images of the modal noise analysis
    radius_factor : float, optional
        fraction of the radius inside which the modal noise is calculated

    Returns
    -------
    signal_to_noise : float
        MEDIAN / STDEV for the intensities inside the fiber face
    """
    image, center, radius = _get_image_data(image_obj, **kwargs)
    if radius_factor is None:
        radius_factor = _get_radius_factor(radius, image_obj.camera)
    inten_array = intensity_array(image, center, radius*radius_factor)

    if show_image:
        tophat_fit = image_obj.get_tophat_fit()
        image, new_center = crop_image(image, center, radius*radius_factor)
        tophat_fit = crop_image(tophat_fit, center, radius*radius_factor)[0]
        plot_overlaid_cross_sections(image, tophat_fit, new_center)
        plot_image(tophat_fit)
        show_plots()

    if image_obj.camera == 'ff':
        return np.max(inten_array) / np.std(inten_array)
    return np.median(inten_array) / np.std(inten_array)

def _modal_noise_contrast(image_obj, radius_factor=None, show_image=False, **kwargs):
    """Finds modal noise of image using Michelson contrast

    Modal noise is defined as (I_max - I_min) / (I_max + I_min)

    Args
    ----
    image_obj : FiberImage
        image object to analyze
    radius_factor : float, optional
        fraction of the radius inside which the modal noise is calculated

    Returns
    -------
    michelson_contrast : float
        (I_max - I_min) / (I_max + I_min) for intensities inside fiber face
    """
    image, center, radius = _get_image_data(image_obj, **kwargs)
    if radius_factor is None:
        if image_obj.camera == 'ff':
            radius_factor = 0.1
        else:
            radius_factor = _get_radius_factor(radius)
    inten_array = intensity_array(image, center, radius*radius_factor)

    return (np.max(inten_array) - np.min(inten_array)) \
           / (np.max(inten_array) + np.min(inten_array))

def _modal_noise_gradient(image_obj, show_image=False, radius_factor=None, **kwargs):
    """Finds modal noise of image using the image gradient

    Modal noise is defined as the stadard deviation of the gradient of the
    fiber face normalized by the median of the original fiber image

    Args
    ----
    image_obj : FiberImage
        image object to analyze
    show_image : bool, optional
        whether or not to show images of the modal noise analysis
    radius_factor : float, optional
        fraction of the radius inside which the modal noise is calculated
    fiber_method : str, optional
        method to use when calculating center and radius of fiber face

    Returns
    -------
    modal_noise : float
        STDEV / MEDIAN for the gradient in the fiber image
    """
    image, center, radius = _get_image_data(image_obj, **kwargs)
    if radius_factor is None:
        radius_factor = _get_radius_factor(radius, image_obj.camera)
    image, center = crop_image(image, center, radius)

    gradient_y, gradient_x = np.gradient(image)
    gradient_array = np.sqrt(gradient_x**2 + gradient_y**2)

    if show_image:
        plot_image(gradient_array)
        plot_overlaid_cross_sections(image, gradient_array, center)
        show_plots()

    inten_array = intensity_array(gradient_array, center, radius*radius_factor)
    image_inten_array = intensity_array(image, center, radius*radius_factor)
    return np.std(inten_array) / np.median(image_inten_array)

def _modal_noise_polynomial(image_obj, show_image=False, radius_factor=None, deg=6, **kwargs):
    """Finds SNR of image using polynomial fit

    Crops image exactly around the circumference of the circle and fits a
    polynomial to the 2D image. Modal noise is then defined as the STDEV
    in the difference between the original image and the polynomial fit
    divided into the median value inside the fiber face

    Args
    ----
    image_obj : FiberImage
        image object to analyze
    show_image : bool, optional
        whether or not to show images of the modal noise analysis
    radius_factor : float, optional
        fraction of the radius inside which the modal noise is calculated
    deg : float, optional
        degree of the fitted polynomial

    Returns
    -------
    signal_to_noise : float
        STDEV for the difference between the fiber image and poly_fit
        divided into the median of the fiber image intensities
    """
    image, center, radius = _get_image_data(image_obj, **kwargs)
    if radius_factor is None:
        radius_factor = _get_radius_factor(radius, image_obj.camera)
    poly_fit = image_obj.get_polynomial_fit(deg, radius_factor, **kwargs)

    if show_image:
        plot_overlaid_cross_sections(image, poly_fit, center)
        plot_image(poly_fit)
        show_plots()

    diff_array = image - poly_fit
    inten_array = intensity_array(diff_array, center, radius * radius_factor)
    poly_inten_array = intensity_array(poly_fit, center, radius * radius_factor)
    # image_inten_array = intensity_array(image, center, radius * radius_factor)
    # return np.median(image_inten_array) / inten_array.std()
    if image_obj.camera == 'ff':
        return np.max(poly_inten_array) / np.std(inten_array)
    return np.median(poly_inten_array) / np.std(inten_array)

def _modal_noise_gaussian(image_obj, show_image=False, radius_factor=None, **kwargs):
    """Finds modal noise of image using a gaussian fit

    Crops image exactly around the circumference of the circle and fits a
    gaussian to the 2D image. Modal noise is then defined as the STDEV
    in the difference between the original image and the gaussian fit
    ndevided into the median value inside the fiber face

    Args
    ----
    image_obj : FiberImage
        image object to analyze
    show_image : bool, optional
        whether or not to show images of the modal noise analysis
    radius_factor : float, optional
        fraction of the radius inside which the modal noise is calculated

    Returns
    -------
    signal_to_noise : float
        STDEV for the difference between the fiber image and gauss_fit
        divided into the median of the fiber image intensities
    """
    image, center, radius = _get_image_data(image_obj, **kwargs)
    if radius_factor is None:
        radius_factor = _get_radius_factor(radius, image_obj.camera)

    gauss_fit = image_obj.get_gaussian_fit()

    if show_image:
        plot_overlaid_cross_sections(image, gauss_fit, center)
        plot_image(gauss_fit)
        show_plots()

    diff_array = image - gauss_fit
    inten_array = intensity_array(diff_array, center, radius*np.sqrt(2))
    gauss_inten_array = intensity_array(gauss_fit, center, radius*np.sqrt(2))
    # image_inten_array = intensity_array(image, center, radius*np.sqrt(2))
    # np.mean(image_inten_array) / return np.std(inten_array)
    return np.max(gauss_inten_array) / np.std(inten_array)

def _modal_noise_gini(image_obj, show_image=False, radius_factor=None, **kwargs):
    """Find modal noise of image using Gini coefficient

    Args
    ----
    image_obj : FiberImage
        image object to analyze
    show_image : bool, optional
        whether or not to show images of the modal noise analysis
    radius_factor : float, optional
        fraction of the radius inside which the modal noise is calculated
    fiber_method : str, optional
        method to use when calculating center and radius of fiber face

    Returns
    -------
    gini_coefficient : float
        Gini coefficient for intensities inside the fiber face
    """
    image, center, radius = _get_image_data(image_obj, **kwargs)
    if radius_factor is None:
        radius_factor = _get_radius_factor(radius, image_obj.camera)

    return _gini_coefficient(intensity_array(image, center,
                                             radius*radius_factor))

def _gini_coefficient(test_array):
    """Finds gini coefficient for intensities in given array

    Args
    ----
    test_array : numpy array
        arbitrary-dimension numpy array of values

    Returns
    -------
    gini_coefficient : float
        Sum( |I_i - I_j| {(i,j), n} ) / ( 2 * n * Sum( I_i, {i, n} ) )
    """
    test_array = test_array.ravel()
    gini_coeff = 0.0
    for i in xrange(test_array.size):
        gini_coeff += np.abs(test_array[i] - test_array).sum()
    norm = 2.0 * test_array.sum() * test_array.size

    return gini_coeff / norm

def _modal_noise_entropy(image_obj, show_image=False, radius_factor=None, **kwargs):
    """Find modal noise of image using hartley entropy

    Args
    ----
    image_obj : FiberImage
        image object to analyze
    show_image : bool, optional
        whether or not to show images of the modal noise analysis
    radius_factor : float, optional
        fraction of the radius inside which the modal noise is calculated

    Returns
    -------
    hartley_entropy : float
        Hartley entropy for intensities inside fiber face
    """
    image, center, radius = _get_image_data(image_obj, **kwargs)
    if radius_factor is None:
        radius_factor = _get_radius_factor(radius, image_obj.camera)

    inten_array = intensity_array(image, center, radius*radius_factor)
    inten_array = inten_array / np.sum(inten_array)
    return np.sum(-inten_array * np.log10(inten_array))

#=============================================================================#
#==== Modal Noise Test =======================================================#
#=============================================================================#

if __name__ == '__main__':
    from fiber_properties import FiberImage, image_list, plot_fft

    base_folder = '../data/Modal Noise Measurements/2016-07-26/'
    ambient_folder = base_folder + 'ambient/600um/'
    dark_folder = base_folder + 'dark/'
    agitated_folder = base_folder + 'images/600um/agitated/'
    unagitated_folder = base_folder + 'images/600um/unagitated/'

    nf = {}
    ff = {}

    nf_dark = image_list(dark_folder + 'nf_dark_')
    nf_ambient = image_list(ambient_folder + 'nf_ambient_')
    ff_dark = image_list(dark_folder + 'ff_dark_')
    ff_ambient = image_list(ambient_folder + 'ff_ambient_')

    nf_agitated = FiberImage(image_list(agitated_folder + 'nf_agitated_'),
                                nf_dark, nf_ambient)
    nf_unagitated = FiberImage(image_list(unagitated_folder + 'nf_unagitated_'),
                                  nf_dark, nf_ambient)
    nf_baseline = FiberImage(nf_agitated.get_tophat_fit(),
                                pixel_size=nf_agitated.get_pixel_size(),
                                threshold=0.1, camera='nf')

    ff_agitated = FiberImage(image_list(agitated_folder + 'ff_agitated_'),
                                ff_dark, ff_ambient)
    ff_unagitated = FiberImage(image_list(unagitated_folder + 'ff_unagitated_'),
                                  ff_dark, ff_ambient)
    ff_baseline = FiberImage(ff_agitated.get_gaussian_fit(),
                                pixel_size=ff_agitated.get_pixel_size(),
                                magnification=1, camera='ff')

    print

    modal_noise = []
    for test in [nf_agitated, nf_unagitated, nf_baseline]:
        modal_noise.append(modal_noise(test, method='fft', output='array', radius_factor=1.0))
    plot_fft([modal_noise[i][1] for i in xrange(3)],
             [modal_noise[i][0] for i in xrange(3)],
             labels=['Agitated laser', 'Unagitated laser', 'Baseline'],
             title='NF Modal Noise Comparison (600um Fiber)')

    modal_noise = []
    for test in [ff_agitated, ff_unagitated, ff_baseline]:
        modal_noise.append(modal_noise(test, method='fft', output='array', radius_factor=1.0))
    plot_fft([modal_noise[i][1] for i in xrange(3)],
             [modal_noise[i][0] for i in xrange(3)],
             labels=['Agitated laser', 'Unagitated laser', 'Baseline'],
             title='FF Modal Noise Comparison (600um Fiber)')

