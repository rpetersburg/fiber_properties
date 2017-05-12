from fiber_properties import (FiberImage, plot_fft, show_plots,
                              save_plot, image_list, filter_image, crop_image)
import numpy as np
import csv
from copy import deepcopy
import time

NEW_DATA = False
NEW_OBJECTS = False
NEW_BASELINE = False
FOLDER = '../data/modal_noise/amp_freq_600um/'
CAMERAS = ['nf', 'ff']
CASE = 3
METHODS = ['tophat', 'gaussian', 'polynomial', 'contrast', 'filter', 'gradient', 'fft']
# METHODS = ['fft']

if CASE == 1:
    TITLE = 'Amplitude vs Frequency'
    TESTS = ['unagitated_10s',
             'agitated_5volts_40mm_10s',
             'agitated_5volts_160mm_10s_test1',
             'agitated_30volts_40mm_10s',
             'agitated_30volts_160mm_10s_test1',
             'baseline_amp_freq']
if CASE == 2:
    TITLE = 'Normalization Test'
    TESTS = ['unagitated_1s',
             'unagitated_8s',
             'unagitated_10s',
             'agitated_5volts_160mm_8s',
             'agitated_5volts_160mm_80s',
             'agitated_30volts_160mm_1s',
             'agitated_30volts_160mm_10s_test2',
             'baseline_norm']
if CASE == 3:
    TITLE = 'Test 1 vs Test 2'
    TESTS = ['unagitated_10s',
             'agitated_5volts_160mm_10s_test1',
             'agitated_5volts_160mm_10s_test2',
             'agitated_30volts_160mm_10s_test1',
             'agitated_30volts_160mm_10s_test2',
             'baseline_norm']

def image_file(test, cam):
    return FOLDER + test + '/' + cam + '_corrected.fit'

def object_file(test, cam):
    return FOLDER + test + '/' + cam + '_obj.pkl'

if __name__ == '__main__':
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
                dark = image_list(FOLDER + 'dark/' + cam + '_')

                ambient_folder = 'ambient_10s/'
                if '1s' in test:
                  ambient_folder = 'ambient_1s/'
                elif '8s' in test:
                  ambient_folder = 'ambient_8s/'
                elif '80s' in test:
                  ambient_folder = 'ambient_80s/'
                ambient = image_list(FOLDER + ambient_folder + cam + '_')
                im_obj = FiberImage(images, dark=dark, ambient=ambient, camera=cam)
                im_obj.save_image(image_file(test, cam))
                im_obj.save_object(object_file(test, cam))

            if NEW_DATA or NEW_OBJECTS:
                im_obj = FiberImage(object_file(test, cam))
                print 'setting new data'
                for method in methods:
                    print 'setting method ' + method
                    im_obj.set_modal_noise(method)
                im_obj.save_object(object_file(test, cam))
            print

        if base_i is not None:
            if NEW_BASELINE:
                print 'saving new baseline object'
                im_obj = FiberImage(object_file(TESTS[base_i-1], cam))
                radius = im_obj.get_fiber_radius()
                kernel_size = int(radius/5.0)
                kernel_size += 1 - (kernel_size % 2)
                image_crop = crop_image(im_obj.get_image(),
                                        im_obj.get_fiber_center(),
                                        radius+kernel_size, False)

                start = time.time()
                perfect_image = filter_image(image_crop, kernel_size=kernel_size)
                perfect_image *= (perfect_image > 0.0).astype('float64')
                end = time.time()
                print('that took ' + str(end-start) + ' seconds')

                baseline_image = np.zeros_like(perfect_image)
                for i in xrange(10):
                    baseline_image += np.random.poisson(perfect_image) / 10.0

                baseline_obj = FiberImage(baseline_image,
                                          pixel_size=im_obj.pixel_size,
                                          camera=cam)
                baseline_obj.save_image(image_file(TESTS[base_i], cam))
                baseline_obj.save_object(object_file(TESTS[base_i], cam))

            if NEW_DATA or NEW_BASELINE:
                baseline_obj = FiberImage(object_file(TESTS[base_i], cam))
                for method in methods:
                    print 'setting method ' + method
                    baseline_obj.set_modal_noise(method)
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

