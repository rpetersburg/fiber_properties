from fiber_properties import (FiberImage, plot_fft, show_plots,
                              save_plot, image_list, filter_image, crop_image)
import numpy as np
import csv
from copy import deepcopy

NEW_DATA = False
NEW_BASELINE = False
CAMERAS = ['ff']
CASE = 5
# METHODS = ['tophat', 'gaussian', 'polynomial', 'contrast', 'filter', 'gradient', 'fft']
METHODS = ['fft']

if CASE == 1:
    TITLE = 'Modal Noise 200-200um'
    FOLDER = '../data/modal_noise/coupled_fibers/200-200um_test2/'
if CASE == 2:
    TITLE = 'Modal Noise 100-200um'
    FOLDER = '../data/modal_noise/coupled_fibers/100-200um/'
if CASE == 3:
    TITLE = 'Modal Noise Octagonal 100um'
    FOLDER = "../data/modal_noise/Kris_data/octagonal_100um/"
if CASE == 4:
    TITLE = 'Modal Noise Circular 100um'
    FOLDER = "../data/modal_noise/Kris_data/circular_100um/"
if CASE == 5:
    TITLE = 'Modal Noise Rectangular 100x300um'
    FOLDER = "../data/modal_noise/Kris_data/rectangular_100x300um/"
if CASE == 6:
    TITLE = 'Modal Noise Unagitated'
    FOLDER = "../data/modal_noise/Kris_data/"
if CASE == 7:
    TITLE = 'Modal Noise Coupled Agitation'
    FOLDER = '../data/modal_noise/Kris_data/'

if CASE == 1:
    TESTS = ['unagitated',
             'agitated_first',
             'agitated_second',
             'agitated_both',
             'baseline']
if CASE == 2:
    TESTS = ['unagitated',
             'agitated_first_100um',
             'agitated_second_200um',
             'agitated_both',
             'baseline']    
if CASE == 3 or CASE == 4 or CASE == 5:
    TESTS = ['unagitated',
             'linear_agitation',
             'circular_agitation',
             'coupled_agitation',
             'baseline']
if CASE == 6:
    TESTS = ['octagonal_100um/unagitated',
             'circular_100um/unagitated',
             'rectangular_100x300um/unagitated']
if CASE == 7:
    TESTS = ['octagonal_100um/coupled_agitation',
             'circular_100um/coupled_agitation',
             'rectangular_100x300um/coupled_agitation']

if CASE == 6 or CASE == 7:
    CAL = ['octagonal_100um/',
           'circular_100um/',
           'rectangular_100x300um']
else:
    CAL = ['' for i in xrange(len(TESTS))]

if __name__ == '__main__':
    for cam in CAMERAS:
        methods = deepcopy(METHODS)
        if cam == 'nf' and 'gaussian' in METHODS:
            methods.remove('gaussian')
        elif cam == 'ff' and 'tophat' in METHODS:
            methods.remove('tophat')

        for i, test in enumerate(TESTS):
            if 'baseline' in test:
                base_i = i
                continue       

            if NEW_DATA:
                print 'saving', cam, test
                images = image_list(FOLDER + test + '/' + cam + '_')
                dark = image_list(FOLDER + CAL[i] + 'dark/' + cam + '_')
                ambient = image_list(FOLDER + CAL[i] + 'ambient/' + cam + '_')
                im_obj = FiberImage(images, dark=dark, ambient=ambient, camera=cam)
                fiber_method = 'edge'
                if 'rectang' in FOLDER or test:
                    fiber_method = 'rectangle'
                for method in methods:
                    print 'setting method ' + method
                    im_obj.set_modal_noise(method, fiber_method=fiber_method)
                im_obj.save_image(FOLDER + test + '/' + cam + '_corrected.fit')
                im_obj.save_object(FOLDER + test + '/' + cam + '_obj.pkl')
                print

        if 'baseline' in TESTS:
            fiber_method = 'edge'
            if 'rectang' in FOLDER:
                fiber_method = 'rectangle'
            if NEW_BASELINE:
                print 'saving', cam, 'baseline'
                im_obj = FiberImage(FOLDER + TESTS[base_i-1] + '/' + cam + '_obj.pkl')
                # perfect_image = im_obj.get_polynomial_fit()
                # perfect_image = im_obj.get_filtered_image(kernel_size=65)
                radius = im_obj.get_fiber_radius(method=fiber_method)
                kernel_size = 101
                # kernel_size = int(radius/3.0)
                # kernel_size += 1 - (kernel_size % 2)
                image_crop = crop_image(im_obj.get_image(),
                                        im_obj.get_fiber_center(method=fiber_method),
                                        radius+kernel_size, False)
                perfect_image = filter_image(image_crop, kernel_size=kernel_size)
                perfect_image *= (perfect_image > 0.0).astype('float64')

                baseline_image = np.zeros_like(perfect_image)
                for i in xrange(10):
                    baseline_image += np.random.poisson(perfect_image) / 10.0

                baseline_obj = FiberImage(baseline_image,
                                          pixel_size=im_obj.pixel_size,
                                          camera=cam)
                baseline_obj.save_image(FOLDER + test + '/' + cam + '_corrected.fit')
                baseline_obj.save_object(FOLDER + test + '/' + cam + '_obj.pkl')

            if NEW_DATA or NEW_BASELINE:
                baseline_obj = FiberImage(FOLDER + 'baseline/' + cam + '_obj.pkl')
                for method in methods:
                    print 'setting method ' + method
                    baseline_obj.set_modal_noise(method, fiber_method=fiber_method)
                baseline_obj.save_image(FOLDER + test + '/' + cam + '_corrected.fit')
                baseline_obj.save_object(FOLDER + test + '/' + cam + '_obj.pkl')

        if 'fft' in methods:
            methods.remove('fft')
            fft_info_list = []
            for test in TESTS:
                im_obj = FiberImage(FOLDER + test + '/' + cam + '_obj.pkl')
                fft_info_list.append(im_obj.get_modal_noise(method='fft'))
            min_wavelength = im_obj.pixel_size / im_obj.magnification * 2.0
            max_wavelength = im_obj.get_fiber_radius(method='edge', units='microns')
            plot_fft(fft_info_list,
                     labels=TESTS,
                     title=TITLE + ' ' + cam.upper(),
                     min_wavelength=min_wavelength,
                     max_wavelength=max_wavelength)
            save_plot(FOLDER + TITLE + ' ' + cam.upper() + '.png')

        modal_noise_info = [['cam', 'test'] + methods]
        for i, test in enumerate(TESTS):
            im_obj = FiberImage(FOLDER + test + '/' + cam + '_obj.pkl')
            modal_noise_info.append([cam, test])         
            for method in methods:
                modal_noise = im_obj.get_modal_noise(method)
                print cam, test, method, modal_noise
                modal_noise_info[i+1].append(modal_noise)

        with open(FOLDER + TITLE + ' ' + cam.upper() + ' Data.csv', 'wb') as f:
            wr = csv.writer(f)
            wr.writerows(modal_noise_info)

