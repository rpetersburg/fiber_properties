from modal_noise_script import (save_new_object, set_new_data,
                                save_baseline_object, save_fft_plot,
                                save_modal_noise_data,
                                save_modal_noise_bar_plot)
from copy import deepcopy
import os

NEW_DATA = False
NEW_OBJECTS = False
NEW_BASELINE = False
OVERWRITE = 'choose'
FOLDER = "C:/Libraries/Box Sync/ExoLab/Fiber_Characterization/Image Analysis/data/modal_noise/rec_fiber_amp_tests/"
CAMERAS = ['nf']
FIBER_METHOD = 'full'
KERNEL = 51
CASE = 1
METHOD = 'filter'
AMBIENT = '../ambient_2s'

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
            new_object = NEW_OBJECTS or cam + '_obj.pkl' not in os.listdir(folder)
            if new_object:
                for i in xrange(10):
                    save_new_object(folder, cam, ambient, num=1, start=i, overwrite=OVERWRITE)
                    save_new_object(folder, cam, ambient, num=i+1, start=0, overwrite=OVERWRITE)

            if NEW_DATA or new_object:
                for i in xrange(10):
                    set_new_data(folder, cam, method, num=1, start=i,
                                 fiber_method=FIBER_METHOD,
                                 kernel_size=KERNEL,
                                 overwrite=OVERWRITE)
                    set_new_data(folder, cam, method, num=i+1, start=0,
                                 fiber_method=FIBER_METHOD,
                                 kernel_size=KERNEL,
                                 overwrite=OVERWRITE)
        if 'fft' in method:
            save_fft_plot(FOLDER, TESTS, cam, LABELS, TITLE)
        else:
            save_modal_noise_bar_plot(FOLDER, TESTS, cam, LABELS, method, TITLE + ' 1x', num=1)
            save_modal_noise_bar_plot(FOLDER, TESTS, cam, LABELS, method, TITLE + ' 10x', num=10)
            save_modal_noise_line_plot(FOLDER, TESTS, cam, LABELS, method, TITLE)
        

if __name__ == '__main__':
    main()
