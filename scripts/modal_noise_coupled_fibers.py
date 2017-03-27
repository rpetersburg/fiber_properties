from fiber_properties import (FiberImage, modal_noise, plot_fft, show_plots,
                              save_plot, image_list)
import numpy as np

FOLDER = '../data/modal_noise/coupled_fibers/2017-03-21/'
CAMERAS = ['nf']
TESTS = ['unagitated',
         'agitated_first',
         'agitated_second',
         'agitated_both',
         'baseline']

if __name__ == '__main__':
    fft = [[] for i in xrange(len(TESTS))]
    freq = [[] for i in xrange(len(TESTS))]
    for cam in CAMERAS:
        for i, test in enumerate(TESTS):
            if 'baseline' in test:
                base_i = i
                continue       

            images = image_list(FOLDER + test + '/' + cam + '_')
            dark = image_list(FOLDER + 'dark/' + cam + '_')
            ambient = image_list(FOLDER + 'ambient/' + cam + '_')
            im_obj = FiberImage(images, dark=dark, ambient=ambient, camera=cam)
            im_obj.save_image(FOLDER + test + '/' + cam + '_corrected.fit')
            
            fft[i], freq[i] = modal_noise(im_obj, method='fft',
                                          output='array', radius_factor=1.0,
                                          show_image=False)

        if 'baseline' in TESTS:
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
            fft[base_i], freq[base_i] = modal_noise(baseline_obj,
                                                    method='fft',
                                                    output='array',
                                                    radius_factor=1.0,
                                                    show_image=False)

        plot_fft(freq, fft,
                 labels=TESTS,
                 title=cam.upper() + ' Modal Noise Coupled Fibers',
                 min_wavelength=1.0,
                 max_wavelength=100.0)
        save_plot(FOLDER + 'Modal Noise.png')
        show_plots()

