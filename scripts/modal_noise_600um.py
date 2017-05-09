from fiber_properties import (FiberImage, plot_fft, show_plots,
                              save_plot, image_list, filter_image, crop_image)
import numpy as np
import csv
from copy import deepcopy

NEW_DATA = False
NEW_BASELINE = False
FOLDER = '../data/modal_noise/amp_freq_600um/'
CAMERAS = ['nf','ff']
CASE = 1
METHODS = ['tophat', 'gaussian', 'polynomial', 'contrast', 'gradient']
# METHODS = ['fft']

if CASE == 1:
    TITLE = 'Amplitude vs Frequency'
    TESTS = ['unagitated_10s',
             'agitated_5volts_40mm_10s',
             'agitated_5volts_160mm_10s_test1',
             'agitated_30volts_40mm_10s',
             'agitated_30volts_160mm_10s_test1']
if CASE == 2:
    TITLE = 'Normalization Test'
    TESTS = ['unagitated_1s',
             'unagitated_8s',
             'unagitated_10s',
             'agitated_5volts_160mm_8s',
             'agitated_5volts_160mm_80s',
             'agitated_30volts_160mm_1s',
             'agitated_30volts_160mm_10s_test2']
if CASE == 3:
    TITLE = 'Test 1 vs Test 2'
    TESTS = ['unagitated_10s',
             'agitated_5volts_160mm_10s_test1',
             'agitated_5volts_160mm_10s_test2',
             'agitated_30volts_160mm_10s_test1',
             'agitated_30volts_160mm_10s_test2']

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
                for method in methods:
                    print 'setting method ' + method
                    im_obj.set_modal_noise(method)
                im_obj.save_image(FOLDER + test + '/' + cam + '_corrected.fit')
                im_obj.save_object(FOLDER + test + '/' + cam + '_obj.pkl')
                print

        if 'baseline' in TESTS:
            if NEW_BASELINE:
                print 'saving', cam, 'baseline'
                im_obj = FiberImage(FOLDER + TESTS[base_i-1] + '/' + cam + '_obj.pkl')
                # perfect_image = im_obj.get_polynomial_fit()
                # perfect_image = im_obj.get_filtered_image(kernel_size=65)
                radius = im_obj.get_fiber_radius()
                kernel_size = int(radius/3.0)
                kernel_size += 1 - (kernel_size % 2)
                image_crop = crop_image(im_obj.get_image(),
                                        im_obj.get_fiber_center(),
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
                    baseline_obj.set_modal_noise(method)
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
                im_obj.save_object(FOLDER + test + '/' + cam + '_obj.pkl')
                print cam, test, method, modal_noise
                modal_noise_info[i+1].append(modal_noise)

        with open(FOLDER + TITLE + ' ' + cam.upper() + ' Data.csv', 'wb') as f:
            wr = csv.writer(f)
            wr.writerows(modal_noise_info)

