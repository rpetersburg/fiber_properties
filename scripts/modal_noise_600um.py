from fiber_properties import (FiberImage, modal_noise, plot_fft, show_plots,
                              save_plot, show_image_array, image_list)
import numpy as np
import re

NEW_DATA = True
FOLDER = '../data/modal_noise/amp_freq_600um/'
CAMERAS = ['nf', 'ff']
# TESTS = ['normalized/unagitated_1s',
#          'normalized/unagitated_8s',
#          'normalized/agitated_5volts_160mm',
#          'normalized/agitated_30volts_160mm',
#          'baseline']
TESTS = ['test1/agitated_5volts_160mm',
         'test1/agitated_30volts_160mm',
         'test2/agitated_5volts_160mm',
         'test2/agitated_30volts_160mm']
TITLE = 'Modal Noise Repeatability Test'

if __name__ == '__main__':
    for cam in CAMERAS:
        if NEW_DATA:
            for i, test in enumerate(TESTS):
                if 'baseline' in test:
                    base_i = i
                    continue

                print 'saving', cam, test
                images = image_list(FOLDER + cam + '/' + test + '_')
                exp_time = int(round(FiberImage(images).exp_time))
                dark = image_list(FOLDER + cam + '/dark/dark_')
                ambient = image_list(FOLDER + cam + '/ambient/ambient_'
                                     + str(exp_time) + 's_')
                im_obj = FiberImage(images, dark=dark, ambient=ambient, camera=cam)
                im_obj.save_image(FOLDER + cam + '/' + test + '_corrected.fit')
                im_obj.set_modal_noise()
                im_obj.save_object(FOLDER + cam + '/' + test + '_obj.pkl')

            if 'baseline' in TESTS:
                    print 'saving', cam, 'baseline'
                    if cam == 'nf':
                        perfect_image = im_obj.get_tophat_fit()
                    elif cam == 'ff':
                        perfect_image = im_obj.get_gaussian_fit()
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
            im_obj = FiberImage(FOLDER + cam + '/' + test + '_obj.pkl')
            fft_info_list.append(im_obj.get_modal_noise(method='fft'))
            if cam == 'nf':
                    method1 = 'tophat'
            elif cam == 'ff':
                method1 = 'gaussian'
            for method in [method1, 'contrast']:
                print cam, test, method, im_obj.get_modal_noise(method)

        min_wavelength = im_obj.pixel_size / im_obj.magnification * 2.0
        max_wavelength = im_obj.get_fiber_radius(method='edge', units='microns')
        plot_fft(fft_info_list,
                 labels=TESTS,
                 title=cam.upper() + TITLE,
                 min_wavelength=min_wavelength,
                 max_wavelength=max_wavelength)
        save_plot(FOLDER + cam + '/' + TITLE + '.png')

