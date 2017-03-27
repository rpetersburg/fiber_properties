from fiber_properties import (FiberImage, modal_noise, plot_fft, show_plots,
                              save_plot, show_image_array, image_list)
import numpy as np
import re

if __name__ == '__main__':
    base_folder = '../data/modal_noise/amp_freq_600um/'
    ext = '.fit'
    cameras = ['nf', 'ff']
    tests = ['normalized/unagitated_1s',
             'normalized/unagitated_8s',
             'normalized/agitated_5volts_160mm',
             'normalized/agitated_30volts_160mm',
             'baseline']

    fft = dict.fromkeys(cameras, {})
    freq = dict.fromkeys(cameras, {})
    labels = []

    for camera in cameras:
        camera_folder = base_folder + camera + '/'
        dark_folder = camera_folder + 'dark/'
        ambient_folder = camera_folder + 'ambient/'

        for test in tests:
            if 'baseline' in test:
                continue

            if '8s' in test:
                suffix = '8s'
            elif '1s' in test:
                suffix = '1s'
            elif 'normalized' in test:
                if '80s' in test:
                    suffix = '80s'
                elif '30volts' in test:
                    suffix = '1s'
                elif '5volts' in test:
                    suffix = '8s'
            else:
                suffix = '10s'

            print test
            images = image_list(camera_folder + test + '_')
            dark = image_list(dark_folder + 'dark_')
            ambient = image_list(ambient_folder + 'ambient_' + suffix + '_')
            im_obj = FiberImage(images, dark, ambient, camera=camera)
            im_obj.save_image()
            
            fft[camera][test], freq[camera][test] = modal_noise(im_obj,
                                                                method='fft',
                                                                output='array',
                                                                radius_factor=1.0,
                                                                show_image=False)
            string_list = re.split('_|/', test)
            if 'agitated' in string_list:
                label = '_'.join(string_list[1:4]) + '_' + str(im_obj.get_image_info('exp_time')) + 's'
            elif 'unagitated' in string_list:
                label = string_list[1] + '_' + str(im_obj.get_image_info('exp_time')) + 's'
            labels.append(label)

        if 'baseline' in tests:
            if camera == 'nf':
                perfect_image = im_obj.get_tophat_fit()
            elif camera == 'ff':
                perfect_image = im_obj.get_gaussian_fit()
                perfect_image *= (perfect_image > 0.0).astype('float64')

            baseline_image = np.zeros_like(perfect_image)
            for i in xrange(10):
                baseline_image += np.random.poisson(perfect_image) / 10

            baseline_obj = FiberImage(baseline_image,
                                         pixel_size=im_obj.get_pixel_size(),
                                         camera=camera)
            fft[camera]['baseline'], freq[camera]['baseline'] = modal_noise(baseline_obj,
                                                                           method='fft',
                                                                           output='array',
                                                                           radius_factor=1.0,
                                                                           show_image=False)
            labels.append('baseline')

        plot_fft([freq[camera][test] for test in tests],
                [fft[camera][test] for test in tests],
                labels=labels,
                title=camera.upper() + ' Modal Noise Comparison (600um Fiber)',
                min_wavelength=1.0,
                max_wavelength=100.0)
        save_plot(camera_folder + 'modal_noise_wavelength')
        show_plots()

