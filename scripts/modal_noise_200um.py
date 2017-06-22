from modal_noise_script import (save_new_object, set_new_data,
                                save_baseline_object, save_fft_plot,
                                save_modal_noise_data)
from copy import deepcopy
import os

NEW_DATA = True
NEW_OBJECTS = True
NEW_BASELINE = False
FOLDER = "C:/Libraries/Box Sync/ExoLab/Fiber_Characterization/Image Analysis/data/modal_noise/amp_freq_200um/"
CAMERAS = ['ff']
FIBER_METHOD = 'edge'
KERNEL = 51
CASE = 9
METHODS = ['tophat', 'gaussian', 'polynomial', 'contrast', 'filter', 'gradient', 'fft']
# METHODS = ['fft']

if CASE == 1:
    TITLE = 'Amplitude vs Frequency'
    TESTS = ['unagitated_1s',
             'agitated_5volts_40mm_1s',
             'agitated_5volts_160mm_1s',
             'agitated_30volts_40mm_1s',
             'agitated_30volts_160mm_1s',
             'baseline']
    LABELS = ['unagitated',
              '0.1Hz 40mm agitation',
              '0.1Hz 160mm agitation',
              '1.0Hz 40mm agitation',
              '1.0Hz 160mm agitation',
              'baseline']
if CASE == 2:
    TITLE = 'Normalization Test'
    TESTS = ['unagitated_8s',
             'unagitated_1s',
             'agitated_5volts_40mm_8s',
             'agitated_5volts_160mm_8s',
             'agitated_30volts_40mm_1s',
             'agitated_30volts_160mm_1s',
             'baseline']
    LABELS = ['unagitated 8s',
              'unagitated 1s',
              '0.1Hz 40mm 8s',
              '0.1Hz 160mm 8s',
              '1.0Hz 40mm 1s',
              '1.0Hz 160mm 1s',
              'LED source']

if CASE == 3:
    TITLE = 'Normalization Test x10'
    TESTS = ['agitated_5volts_40mm_80s',
             'agitated_5volts_160mm_80s',
             'agitated_30volts_40mm_10s',
             'agitated_30volts_160mm_10s']
    LABELS = ['0.1Hz 40mm 80s',
              '0.1Hz 160mm 80s',
              '1.0Hz 40mm 10s',
              '1.0Hz 160mm 10s']

if CASE == 4:
    TITLE = 'Exposure Test'
    TESTS = ['agitated_5volts_40mm_1s',
             'agitated_5volts_40mm_8s',
             'agitated_5volts_40mm_80s',
             'agitated_5volts_160mm_1s',
             'agitated_5volts_160mm_8s',
             'agitated_5volts_160mm_80s']
    LABELS = ['0.1Hz 40mm 1s',
              '0.1Hz 40mm 8s',
              '0.1Hz 40mm 80s',
              '0.1Hz 160mm 1s',
              '0.1Hz 160mm 8s',
              '0.1Hz 160mm 80s']

if CASE == 5:
    TITLE = 'Exposure Test 2'
    TESTS = ['agitated_30volts_40mm_1s',
             'agitated_30volts_40mm_10s',
             'agitated_30volts_160mm_1s',
             'agitated_30volts_160mm_10s']
    LABELS = ['1.0Hz 40mm 1s',
              '1.0Hz 40mm 10s',
              '1.0Hz 160mm 1s',
              '1.0Hz 160mm 10s']

if CASE == 6:
    TITLE = 'Exposure Test 3'
    TESTS = ['agitated_5volts_160mm_1s',
             'agitated_5volts_160mm_8s',
             'agitated_5volts_160mm_80s',
             'agitated_30volts_160mm_1s',
             'agitated_30volts_160mm_10s']
    LABELS = ['0.1Hz 160mm 1s',
              '0.1Hz 160mm 8s',
              '0.1Hz 160mm 80s',
              '1.0Hz 160mm 1s',
              '1.0Hz 160mm 10s']

if CASE == 7:
    TITLE = 'Exposure Test 4'
    TESTS = ['agitated_5volts_40mm_1s',
             'agitated_5volts_40mm_8s',
             'agitated_5volts_40mm_80s',
             'agitated_30volts_40mm_1s',
             'agitated_30volts_40mm_10s']
    LABELS = ['0.1Hz 40mm 1s',
              '0.1Hz 40mm 8s',
              '0.1Hz 40mm 80s',
              '1.0Hz 40mm 1s',
              '1.0Hz 40mm 10s']

if CASE == 8:
    TITLE = 'Exposure Test 5'
    TESTS = ['agitated_5volts_40mm_1s',
             'agitated_5volts_40mm_8s',
             'agitated_5volts_40mm_80s',
             'agitated_30volts_40mm_1s',
             'agitated_30volts_40mm_10s',
             'agitated_5volts_160mm_1s',
             'agitated_5volts_160mm_8s',
             'agitated_5volts_160mm_80s',
             'agitated_30volts_160mm_1s',
             'agitated_30volts_160mm_10s',
             'baseline']
    LABELS = ['0.1Hz 40mm 1s',
              '0.1Hz 40mm 8s',
              '0.1Hz 40mm 80s',
              '1.0Hz 40mm 1s',
              '1.0Hz 40mm 10s',
              '0.1Hz 160mm 1s',
              '0.1Hz 160mm 8s',
              '0.1Hz 160mm 80s',
              '1.0Hz 160mm 1s',
              '1.0Hz 160mm 10s',
              'baseline']

if CASE == 9:
    TITLE = 'All'
    TESTS = ['unagitated_1s',
             'unagitated_8s',
             'agitated_5volts_40mm_1s',
             'agitated_5volts_40mm_8s',
             'agitated_5volts_40mm_80s',
             'agitated_30volts_40mm_1s',
             'agitated_30volts_40mm_10s',
             'agitated_5volts_160mm_1s',
             'agitated_5volts_160mm_8s',
             'agitated_5volts_160mm_80s',
             'agitated_30volts_160mm_1s',
             'agitated_30volts_160mm_10s',
             'baseline']
    LABELS = ['0.1Hz 40mm 1s',
              '0.1Hz 40mm 8s',
              '0.1Hz 40mm 80s',
              '1.0Hz 40mm 1s',
              '1.0Hz 40mm 10s',
              '0.1Hz 160mm 1s',
              '0.1Hz 160mm 8s',
              '0.1Hz 160mm 80s',
              '1.0Hz 160mm 1s',
              '1.0Hz 160mm 10s',
              'baseline']

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
            if 'baseline_filter' in test:
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
                set_new_data(FOLDER, test, cam, methods, FIBER_METHOD, KERNEL)

        if base_i is not None:
            new_baseline = NEW_BASELINE or cam + '_obj.pkl' not in os.listdir(FOLDER + TESTS[base_i] + '/')
            if new_baseline:
                save_baseline_object(FOLDER, TESTS[base_i], cam, TESTS[base_i-1], FIBER_METHOD)

            if NEW_DATA or new_baseline:
                set_new_data(FOLDER, TESTS[base_i], cam, FIBER_METHOD, KERNEL)
        if 'fft' in methods:
            methods.remove('fft')
            save_fft_plot(FOLDER, TESTS, cam, LABELS, TITLE)

        save_modal_noise_data(FOLDER, TESTS, cam, methods, TITLE)
