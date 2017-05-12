from fiber_properties import (FiberImage, plot_fft, show_plots,
                              save_plot, image_list, baseline_image)
import numpy as np
import csv
from copy import deepcopy

NEW_DATA = False
NEW_OBJECTS = False
NEW_BASELINE = True
FIBER_METHOD = 'edge'
CAMERAS = ['nf']
CASE = 4
METHODS = ['tophat', 'gaussian', 'polynomial', 'contrast', 'filter', 'gradient', 'fft']
# METHODS = ['filter', 'fft']
# METHODS = ['fft']

if CASE == 1:
    TITLE = 'Modal Noise 200-200um'
    FOLDER = '../data/modal_noise/coupled_fibers/200-200um_test2/'
    TESTS = ['unagitated',
             'agitated_first',
             'agitated_second',
             'agitated_both',
             'baseline']
if CASE == 2:
    TITLE = 'Modal Noise 100-200um'
    FOLDER = '../data/modal_noise/coupled_fibers/100-200um/'
    TESTS = ['unagitated',
             'agitated_first_100um',
             'agitated_second_200um',
             'agitated_both',
             'baseline']
if CASE == 3:
    TITLE = 'Modal Noise Octagonal 100um'
    FOLDER = "../data/modal_noise/Kris_data/octagonal_100um/"
if CASE == 4:
    TITLE = 'Modal Noise Circular 100um'
    FOLDER = "../data/modal_noise/Kris_data/circular_100um/"
if CASE == 5:
    TITLE = 'Modal Noise Rectangular 100x300um'
    FOLDER = "../data/modal_noise/Kris_data/rectangular_100x300um/"
if CASE == 3 or CASE == 4 or CASE == 5:
    TESTS = ['unagitated',
             'linear_agitation',
             'circular_agitation',
             'coupled_agitation',
             'baseline']
if CASE == 6:
    TITLE = 'Modal Noise Unagitated'
    FOLDER = "../data/modal_noise/Kris_data/"
    TESTS = ['octagonal_100um/unagitated',
             'circular_100um/unagitated',
             'rectangular_100x300um/unagitated']
if CASE == 7:
    TITLE = 'Modal Noise Coupled Agitation'
    FOLDER = '../data/modal_noise/Kris_data/'
    TESTS = ['octagonal_100um/coupled_agitation',
             'circular_100um/coupled_agitation',
             'rectangular_100x300um/coupled_agitation']
if CASE == 8:
    TITLE = 'Modal Noise Geometry Baselines'
    FOLDER = '../data/modal_noise/Kris_data/'
    TESTS = ['octagonal_100um/baseline',
             'circular_100um/baseline',
             'rectangular_100x300um/baseline']

if CASE == 6 or CASE == 7:
    CAL = ['octagonal_100um/',
           'circular_100um/',
           'rectangular_100x300um']
else:
    CAL = ['' for i in xrange(len(TESTS))]

def image_file(test, cam):
    return FOLDER + test + '/' + cam + '_corrected.fit'

def object_file(test, cam):
    return FOLDER + test + '/' + cam + '_obj.pkl'

if __name__ == '__main__':
    print TITLE
    print
    for cam in CAMERAS:
        methods = deepcopy(METHODS)
        if cam == 'nf' and 'gaussian' in METHODS:
            methods.remove('gaussian')
        elif cam == 'ff' and 'tophat' in METHODS:
            methods.remove('tophat')

        base_i = None
        for i, test in enumerate(TESTS):
            if 'baseline' in test:
                base_i = i
                continue       

            print cam, test
            if NEW_OBJECTS:
                print 'saving new object'
                images = image_list(FOLDER + test + '/' + cam + '_')
                dark = image_list(FOLDER + CAL[i] + 'dark/' + cam + '_')
                ambient = image_list(FOLDER + CAL[i] + 'ambient/' + cam + '_')
                im_obj = FiberImage(images, dark=dark, ambient=ambient, camera=cam)
                im_obj.save_image(image_file(test, cam))
                im_obj.save_object(object_file(test, cam))

            if NEW_DATA or NEW_OBJECTS:
                print 'setting new data'
                im_obj = FiberImage(object_file(test, cam))
                for method in methods:
                    print 'setting method ' + method
                    im_obj.set_modal_noise(method, fiber_method=FIBER_METHOD)
                im_obj.save_object(object_file(test, cam))
                print

        if base_i is not None:
            if NEW_BASELINE:
                print 'saving new baseline object'
                im_obj = FiberImage(object_file(TESTS[base_i-1], cam))
                baseline_image = baseline_image(im_obj, fiber_method=FIBER_METHOD)

                baseline_obj = FiberImage(baseline_image,
                                          pixel_size=im_obj.pixel_size,
                                          camera=cam)
                baseline_obj.save_image(image_file(TESTS[base_i], cam))
                baseline_obj.save_object(object_file(TESTS[base_i], cam))

            if NEW_DATA or NEW_BASELINE:
                print 'setting new baseline data'
                baseline_obj = FiberImage(FOLDER + 'baseline/' + cam + '_obj.pkl')
                for method in methods:
                    print 'setting method ' + method
                    baseline_obj.set_modal_noise(method, fiber_method=FIBER_METHOD)
                baseline_obj.save_object(object_file(TESTS[base_i], cam))
                print

        if 'fft' in methods:
            print 'saving fft plot'
            methods.remove('fft')
            fft_info_list = []
            for test in TESTS:
                im_obj = FiberImage(FOLDER + test + '/' + cam + '_obj.pkl')
                fft_info_list.append(im_obj.get_modal_noise(method='fft'))
            min_wavelength = im_obj.pixel_size / im_obj.magnification * 2.0
            max_wavelength = im_obj.get_fiber_radius(method='edge', units='microns')
            plot_fft(fft_info_list,
                     labels=TESTS,
                     min_wavelength=min_wavelength,
                     max_wavelength=max_wavelength)
            save_plot(FOLDER + TITLE + ' ' + cam.upper() + '.png')

        print 'saving modal noise data'
        modal_noise_info = [['cam', 'test'] + methods]
        for i, test in enumerate(TESTS):
            im_obj = FiberImage(object_file(test, cam))
            modal_noise_info.append([cam, test])         
            for method in methods:
                modal_noise = im_obj.get_modal_noise(method)
                im_obj.save_object(object_file(test, cam))
                print cam, test, method, modal_noise
                modal_noise_info[i+1].append(modal_noise)

        with open(FOLDER + TITLE + ' ' + cam.upper() + ' Data.csv', 'wb') as f:
            wr = csv.writer(f)
            wr.writerows(modal_noise_info)

