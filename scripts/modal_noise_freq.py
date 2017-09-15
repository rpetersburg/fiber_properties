from modal_noise_script import (save_new_object, set_new_data,
                                save_baseline_object, save_fft_plot,
                                save_modal_noise_data,
                                save_modal_noise_bar_plot,
                                save_modal_noise_inside,
                                save_modal_noise_line_plot)
from copy import deepcopy
import os

# NEW_DATA = False
# NEW_OBJECTS = False
# NEW_BASELINE = False
CAMERAS = ['nf']
FIBER_METHOD = 'full'
KERNEL = 51
CASE = 1
METHOD = 'filter'
OVERWRITE = False
FOLDER = '../data/modal_noise/rec_fiber_freq_tests/'

if CASE == 1:
    FOLDER += "circular_ag/"
    TITLE = 'Circular Frequency'
    TESTS = ['agitated_5volts_6.3s',
             'agitated_10volts_2.6s',
             'agitated_15volts_1.7s',
             'agitated_20volts_1.2s',
             '../unagitated_1s']
    AMBIENTS = ['../ambient_6.3s/',
                '../ambient_2.6s/',
                '../ambient_1.7s/',
                '../ambient_1.2s/',
                '../ambient_1s/']
    LABELS = ['0.16Hz agitation',
              '0.38Hz agitation',
              '0.59Hz agitation',
              '0.83Hz agitation',
              'unagitated']

if CASE == 2:
    TITLE = 'Linear Frequency'
    TESTS = ['agitated_5volts_120mm_7.5s',
             'agitated_10volts_120mm_3.2s',
             'agitated_15volts_120mm_2s',
             'agitated_20volts_120mm_1.5s',
             'agitated_25volts_120mm_1.2s',
             'agitated_30volts_120mm_1s',
             'unagitated_1s']
    AMBIENTS = ['../ambient_7.5s/',
                '../ambient_3.2s',
                '../ambient_2s/',
                '../ambient_1.5s/',
                '../ambient_1.2s/',
                '../ambient_1s/',
                '../ambient_1s/']
    LABELS = ['0.13Hz agitation',
              '0.31Hz agitation',
              '0.50Hz agitation',
              '0.67Hz agitation',
              '0.83Hz agitation',
              '1.00Hz agitation',
              'unagitated']

def main():
    print TITLE
    print
    for cam in CAMERAS:

        kernel = KERNEL
        method = METHOD

        for test, ambient in zip(TESTS, AMBIENTS):
            folder = FOLDER + test
            
            print cam, test
            save_modal_noise_inside(folder, [cam], [method], overwrite=OVERWRITE,
                                    ambient_folder=ambient, kernel_size=KERNEL,
                                    fiber_method=FIBER_METHOD)
            # new_object = NEW_OBJECTS or cam + '_obj.pkl' not in os.listdir(folder)
            # if new_object:
            #     for i in xrange(10):
            #         save_new_object(folder, cam, ambient, num=1, start=i, overwrite=OVERWRITE)
            #         save_new_object(folder, cam, ambient, num=i+1, start=0, overwrite=OVERWRITE)

            # if NEW_DATA or new_object:
            #     for i in xrange(10):
            #         set_new_data(folder, cam, method, num=1, start=i,
            #                      fiber_method=FIBER_METHOD,
            #                      kernel_size=KERNEL,
            #                      overwrite=OVERWRITE)
            #         set_new_data(folder, cam, method, num=i+1, start=0,
            #                      fiber_method=FIBER_METHOD,
            #                      kernel_size=KERNEL,
            #                      overwrite=OVERWRITE)
        if 'fft' in method:
            save_fft_plot(FOLDER, TESTS, cam, LABELS, TITLE)
        else:
            save_modal_noise_bar_plot(FOLDER, TESTS, cam, LABELS, method, TITLE + ' 1x', num=1)
            save_modal_noise_bar_plot(FOLDER, TESTS, cam, LABELS, method, TITLE + ' 10x', num=10)
            save_modal_noise_line_plot(FOLDER, TESTS, cam, LABELS, method, TITLE)
        

if __name__ == '__main__':
    main()
