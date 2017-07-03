from modal_noise_script import (save_new_object, set_new_data,
                                save_baseline_object, save_fft_plot,
                                save_modal_noise_data)
import numpy as np
import os
import csv
from copy import deepcopy

NEW_DATA = True
NEW_OBJECTS = False
NEW_BASELINE = False
FOLDER = "C:/Libraries/Box Sync/ExoLab/Fiber_Characterization/Image Analysis/data/modal_noise/amp_freq_600um/"
CAMERAS = ['ff']
KERNEL = 101
FIBER_METHOD = 'edge'
CASE = 1
# METHODS = ['tophat', 'gaussian', 'polynomial', 'contrast', 'filter', 'gradient', 'fft']
# METHODS = ['tophat', 'gaussian', 'polynomial', 'contrast', 'filter', 'gradient']
METHODS = ['filter', 'fft']

if CASE == 1:
    TITLE = 'Amplitude vs Frequency'
    TESTS = ['unagitated_10s',
             'agitated_5volts_40mm_10s',
             'agitated_5volts_160mm_10s_test1',
             'agitated_30volts_40mm_10s',
             'agitated_30volts_160mm_10s_test1']
    LABELS = ['unagitated',
              '0.1Hz 40mm agitation',
              '0.1Hz 160mm agitation',
              '1.0Hz 40mm agitation',
              '1.0Hz 160mm agitation']
if CASE == 2:
    TITLE = 'Normalization Test'
    TESTS = ['unagitated_1s',
             'unagitated_8s',
             'unagitated_10s',
             'agitated_5volts_160mm_8s',
             'agitated_5volts_160mm_80s',
             'agitated_30volts_160mm_1s',
             'agitated_30volts_160mm_10s_test2']
    LABELS = ['unagitated 1s-exp',
              'unagitated 8s-exp',
              'unagitated 10s-exp',
              '0.1Hz agitation 8s-exp',
              '0.1Hz agitation 80s-exp',
              '1.0Hz agitation 1s-exp',
              '1.0Hz agitation 10s-exp']
if CASE == 3:
    TITLE = 'Test 1 vs Test 2'
    TESTS = ['unagitated_10s',
             'agitated_5volts_160mm_10s_test1',
             'agitated_5volts_160mm_10s_test2',
             'agitated_30volts_160mm_10s_test1',
             'agitated_30volts_160mm_10s_test2']
    LABELS = ['unagitated',
              '0.1Hz agitation test 1',
              '0.1Hz agitation test 2',
              '1.0Hz agitation test 1',
              '1.0Hz agitation test 2']

if CASE == 4:
    TITLE = 'All'
    TESTS = ['unagitated_1s',
             'unagitated_8s',
             'unagitated_10s',
             'agitated_5volts_40mm_10s',
             'agitated_5volts_160mm_8s',
             'agitated_5volts_160mm_10s_test1',
             'agitated_5volts_160mm_10s_test2',
             'agitated_5volts_160mm_80s',
             'agitated_30volts_40mm_10s',
             'agitated_30volts_160mm_1s',
             'agitated_30volts_160mm_10s_test1',
             'agitated_30volts_160mm_10s_test2']
    LABELS = ['unagitated 1s',
              'unagitated 8s',
              'unagitated 10s',
              '0.1Hz 40mm 10s',
              '0.1Hz 160mm 8s',
              '0.1Hz 160mm 10s test1',
              '0.1Hz 160mm 10s test2',
              '0.1Hz 160mm 80s',
              '1.0Hz 40mm 10s',
              '1.0Hz 160mm 1s',
              '1.0Hz 160mm 10s test1',
              '1.0Hz 160mm 10s test2']

if __name__ == '__main__':
    print TITLE
    print
    for cam in CAMERAS:
        methods = deepcopy(METHODS)
        if cam == 'nf' and 'gaussian' in METHODS:
            methods.remove('gaussian')
        elif cam == 'ff' and 'tophat' in METHODS:
            methods.remove('tophat')   

        kernel = KERNEL
        if cam == 'ff':
            kernel = None     

        base_i = None
        for i, test in enumerate(TESTS):
            if 'baseline' in test:
                base_i = i
                continue

            print cam, test
            new_object = NEW_OBJECTS or cam + '_obj.pkl' not in os.listdir(FOLDER + test + '/')
            if new_object:
                dark_folder = 'dark/'
                ambient_folder = 'ambient_1s/'
                if '8s' in test:
                    ambient_folder = 'ambient_8s/'
                if '10s' in test:
                    ambient_folder = 'ambient_10s/'
                if '80s' in test:
                    ambient_folder = 'ambient_80s/'
                save_new_object(FOLDER, test, cam, ambient_folder)

            if NEW_DATA or new_object:
                set_new_data(FOLDER + test, cam, methods,
                             fiber_method=FIBER_METHOD,
                             kernel_size=KERNEL)

        if base_i is not None:
            new_baseline = NEW_BASELINE or cam + '_obj.pkl' not in os.listdir(FOLDER + TESTS[base_i] + '/')
            if new_baseline:
                save_baseline_object(FOLDER + TESTS[base_i], cam,
                                     TESTS[base_i-1],
                                     fiber_method=FIBER_METHOD,
                                     kernel=KERNEL)

            if NEW_DATA or new_baseline:
                set_new_data(FOLDER + TESTS[base_i], cam, methods,
                             fiber_method=FIBER_METHOD,
                             kernel_size=KERNEL)

        if 'fft' in methods:
            methods.remove('fft')
            save_fft_plot(FOLDER, test, cam, LABELS, TITLE)

        save_modal_noise_data(FOLDER, TESTS, cam, LABELS, methods, TITLE)
