from modal_noise_script import (save_fft_plot,
                                save_modal_noise_bar_plot,
                                save_modal_noise_line_plot,
                                save_modal_noise_inside)
from copy import deepcopy
import os

OVERWRITE = False
FOLDER = "../data/modal_noise/rec_fiber_amp_tests/"
CAMERAS = ['nf']
FIBER_METHOD = 'full'
KERNEL = 51
CASE = 1
METHOD = 'filter'

if CASE == 1:
    TITLE = 'Amplitude'
    TESTS = ['agitated_15volts_40mm_2s',
             'agitated_15volts_80mm_2s',
             'agitated_15volts_120mm_2s',
             'agitated_15volts_160mm_2s',
             '../rec_fiber_freq_tests/unagitated_1s']
    AMBIENTS = ['../ambient_2s' for i in xrange(4)]
    AMBIENTS += ['../rec_fiber_freq_tests/ambient_1s']
    LABELS = ['80mm agitation',
              '160mm agitation',
              '240mm agitation',
              '320mm agitation',
              'unagitated']

def main():
    print TITLE
    print

    method = METHOD
    for cam in CAMERAS:

        for test, ambient in zip(TESTS, AMBIENTS):
            folder = FOLDER + test

            print cam, test
            save_modal_noise_inside(folder, [cam], [method], overwrite=OVERWRITE,
                                    ambient_folder=ambient, kernel_size=KERNEL,
                                    fiber_method=FIBER_METHOD)

        if 'fft' in method:
            save_fft_plot(FOLDER, TESTS, cam, LABELS, TITLE)
        else:
            save_modal_noise_bar_plot(FOLDER, TESTS, cam, LABELS, method, TITLE + ' 1x', num=1)
            save_modal_noise_bar_plot(FOLDER, TESTS, cam, LABELS, method, TITLE + ' 10x', num=10)
            save_modal_noise_line_plot(FOLDER, TESTS, cam, LABELS, method, TITLE)
        

if __name__ == '__main__':
    main()
