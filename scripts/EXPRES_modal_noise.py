from fiber_properties import (ImageAnalysis, modal_noise, plot_fft, show_plots,
                              save_plot, show_image_array, image_list)
import numpy as np
import re

FOLDER = '../data/EXPRES/rectangular_132/modal_noise/'
CAMS = ['nf', 'ff']
METHOD = 'rectangle'

if __name__ == '__main__':
    for cam in CAMS:
        images = image_list(FOLDER + cam + '_')
        dark = image_list(FOLDER + '../dark/' + cam + '_')
        ambient = image_list(FOLDER + '../dark/' + cam + '_')
        im_obj = ImageAnalysis(images, dark, ambient, camera=cam)
        im_obj.save_image()

        fft, freq = modal_noise(im_obj,
                                method='fft',
                                output='array',
                                radius_factor=1.0,
                                show_image=False,
                                fiber_method=METHOD)

        plot_fft([freq], [fft],
                 title=cam.upper() + ' Modal Noise',
                 min_wavelength=1.0,
                 max_wavelength=20.0)
        save_plot(FOLDER + cam + '_modal_noise.png')
        show_plots()

