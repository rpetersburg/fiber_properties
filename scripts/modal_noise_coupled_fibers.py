from fiber_properties import (FiberImage, plot_fft, show_plots,
                              save_plot, image_list, filter_image)
import numpy as np
from tabulate import tabulate

NEW_DATA = True
CAMERAS = ['nf']
CASE = 4
TITLE = 'Modal Noise Circular 100um'

if CASE == 1:
    FOLDER = '../data/modal_noise/coupled_fibers/200-200um_test2/'
if CASE == 2:
    FOLDER = '../data/modal_noise/coupled_fibers/100-200um/'
if CASE == 3:
    FOLDER = "../data/modal_noise/Kris' data/octagonal_100um/"
if CASE == 4:
    FOLDER = "../data/modal_noise/Kris' data/circular_100um/"
if CASE == 5:
    FOLDER = "../data/modal_noise/Kris' data/rectangular_100x300um/"
if CASE == 6:
    FOLDER = "../data/modal_noise/Kris' data/"

if CASE == 1 or CASE == 2:
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
             'octagonal_100um/coupled_agitation',
             'circular_100um/unagitated',
             'circular_100um/coupled_agitation']

if CASE == 6:
    CAL = ['octagonal_100um/',
           'octagonal_100um/',
           'circular_100um/',
           'circular_100um/']
else:
    CAL = ['' for i in xrange(len(TESTS))]

if __name__ == '__main__':
    for cam in CAMERAS:
        if NEW_DATA:
            for i, test in enumerate(TESTS):
                if 'baseline' in test:
                    base_i = i
                    continue       

                print 'saving', cam, test
                images = image_list(FOLDER + test + '/' + cam + '_')
                dark = image_list(FOLDER + CAL[i] + 'dark/' + cam + '_')
                ambient = image_list(FOLDER + CAL[i] + 'ambient/' + cam + '_')
                im_obj = FiberImage(images, dark=dark, ambient=ambient, camera=cam)
                im_obj.set_modal_noise(method='filter', show_image=True)
                im_obj.save_image(FOLDER + test + '/' + cam + '_corrected.fit')
                im_obj.save_object(FOLDER + test + '/' + cam + '_obj.pkl')                

            if 'baseline' in TESTS:
                print 'saving', cam, 'baseline'
                # if cam == 'nf':
                #     perfect_image = im_obj.get_tophat_fit()
                # elif cam == 'ff':
                #     perfect_image = im_obj.get_gaussian_fit()
                #     perfect_image *= (perfect_image > 0.0).astype('float64')
                perfect_image = im_obj.get_filtered_image(kernel_size=65)
                perfect_image *= (perfect_image > 0.0).astype('float64')

                baseline_image = np.zeros_like(perfect_image)
                for i in xrange(10):
                    baseline_image += np.random.poisson(perfect_image) / 10.0

                baseline_obj = FiberImage(baseline_image,
                                          pixel_size=im_obj.get_pixel_size(),
                                          camera=cam)
                baseline_obj.set_modal_noise()
                baseline_obj.save_image(FOLDER + test + '/' + cam + '_corrected.fit')
                baseline_obj.save_object(FOLDER + test + '/' + cam + '_obj.pkl')

        fft_info_list = []
        for test in TESTS:
            im_obj = FiberImage(FOLDER + test + '/' + cam + '_obj.pkl')
            fft_info_list.append(im_obj.get_modal_noise(method='fft'))
            if cam == 'nf':
                    method1 = 'tophat'
            elif cam == 'ff':
                method1 = 'polynomial'
            for method in [method1, 'filter', 'contrast']:
                print cam, test, method, im_obj.get_modal_noise(method)

        min_wavelength = im_obj.pixel_size / im_obj.magnification * 2.0
        max_wavelength = im_obj.get_fiber_radius(method='edge', units='microns')
        plot_fft(fft_info_list,
                 labels=TESTS,
                 title=TITLE + ' ' + cam.upper(),
                 min_wavelength=min_wavelength,
                 max_wavelength=max_wavelength)
        save_plot(FOLDER + TITLE + ' ' + cam.upper() + '.png')
        show_plots()

