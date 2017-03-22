from fiber_properties import (ImageAnalysis, modal_noise, plot_fft,
                              show_plots, save_plot, image_list)
import numpy as np

FOLDER = "../data/modal_noise/Kris' data/Circular 100um/"
CAMERAS = ['nf', 'ff']
TESTS = ['unagitated',
         'linear_agitation',
         'circular_agitation',
         'coupled_agitation',
         'baseline']

if __name__ == '__main__':
    fft = [[] for i in xrange(len(TESTS))]
    freq = [[] for i in xrange(len(TESTS))]
    for cam in CAMERAS:
        for i, test in enumerate(TESTS):
            if 'baseline' in test:
                base_i = i
                continue       

            print cam, test
            images = image_list(FOLDER + cam + '/circular_100um_' + test + '_')
            dark = image_list(FOLDER + cam + '/Dark_2_17_2017_')
            ambient = image_list(FOLDER + cam + '/Ambient_2_17_2017_')
            im_obj = ImageAnalysis(images, dark=dark, ambient=ambient, camera=cam)
            im_obj.save_image(FOLDER + cam + '/circular_100um_' + test + '_corrected.fit')
            
            fft[i], freq[i] = modal_noise(im_obj, method='fft',
                                          output='array', radius_factor=1.0,
                                          show_image=False)

        if 'baseline' in TESTS:
            print cam, 'baseline'
            if cam == 'nf':
                perfect_image = im_obj.get_tophat_fit()
            elif cam == 'ff':
                perfect_image = im_obj.get_gaussian_fit()
                perfect_image *= (perfect_image > 0.0).astype('float64')

            baseline_image = np.zeros_like(perfect_image)
            for i in xrange(10):
                baseline_image += np.random.poisson(perfect_image) / 10.0

            baseline_obj = ImageAnalysis(baseline_image,
                                         pixel_size=im_obj.get_pixel_size(),
                                         camera=cam)
            fft[base_i], freq[base_i] = modal_noise(baseline_obj,
                                                    method='fft',
                                                    output='array',
                                                    radius_factor=1.0,
                                                    show_image=False)

        if cam == 'ff':
            min_wavelength = 20.0
            max_wavelength = 5000.0
        else:
            min_wavelength = 1.0
            max_wavelength = 100.0
        plot_fft(freq, fft,
                 labels=TESTS,
                 title=cam.upper() + ' Modal Noise Agitation Method',
                 min_wavelength=min_wavelength,
                 max_wavelength=max_wavelength)
        save_plot(FOLDER + cam.upper() + ' Modal Noise.png')
        show_plots()

