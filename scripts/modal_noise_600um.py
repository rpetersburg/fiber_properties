from modal_noise_script import (save_new_object, set_new_data,
                                save_baseline_object, save_fft_plot,
                                save_modal_noise_data)
import numpy as np
import csv
from copy import deepcopy

NEW_DATA = False
NEW_OBJECTS = False
NEW_BASELINE = False
FOLDER = "C:/Libraries/Box Sync/ExoLab/Fiber_Characterization/Image Analysis/data/modal_noise/amp_freq_600um/"
CAMERAS = ['nf', 'ff']
FIBER_METHOD = 'edge'
CASE = 3
METHODS = ['tophat', 'gaussian', 'polynomial', 'contrast', 'filter', 'gradient', 'fft']
# METHODS = ['fft']

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

if __name__ == '__main__':
    print TITLE
    print
    for cam in CAMERAS:
        methods = deepcopy(METHODS)
        if cam == 'nf' and 'gaussian' in METHODS:
            methods.remove('gaussian')
        elif cam == 'ff' and 'tophat' in METHODS:
            methods.remove('tophat')        

        base_i = None
        for i, test in enumerate(TESTS):
            if 'baseline' in test:
                base_i = i
                continue       

            print cam, test
            new_object = NEW_OBJECTS or cam + '_obj.pkl' not in os.listdir(FOLDER + test + '/')
            if new_object:
                ambient_folder = 'ambient_1s/'
                if '8s' in test:
                    ambient_folder = 'ambient_8s/'
                if '10s' in test:
                    ambient_folder = 'ambient_10s/'
                if '80s' in test:
                    ambient_folder = 'ambient_80s/'
                save_new_object(FOLDER, test, cam, ambient_folder)

            if NEW_DATA or new_object:
                set_new_data(FOLDER, test, cam, methods, FIBER_METHOD)

        if base_i is not None:
            new_baseline = NEW_BASELINE or cam + '_obj.pkl' not in os.listdir(FOLDER + TESTS[base_i] + '/')
            if new_baseline:
                save_baseline_object(FOLDER, TESTS[base_i], cam, TESTS[base_i-1], FIBER_METHOD)

            if NEW_DATA or new_baseline:
                set_new_data(FOLDER, TESTS[base_i], cam, methods, FIBER_METHOD)

        if 'fft' in methods:
            methods.remove('fft')
            save_fft_plot(FOLDER, test, cam, LABELS, TITLE)

        save_modal_noise_data(FOLDER, test, cam, methods, TITLE)
