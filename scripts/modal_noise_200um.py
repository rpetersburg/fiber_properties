from modal_noise_script import (save_new_object, set_new_data,
                                save_baseline_object, save_fft_plot,
                                save_modal_noise_data,
                                save_modal_noise_bar_plot)
from copy import deepcopy
import os

NEW_DATA = False
NEW_OBJECTS = False
NEW_BASELINE = False
FOLDER = "C:/Libraries/Box Sync/ExoLab/Fiber_Characterization/Image Analysis/data/modal_noise/amp_freq_200um/"
CAMERAS = ['nf', 'ff']
FIBER_METHOD = 'edge'
KERNEL = 51
CASE = 8
# METHODS = ['tophat', 'gaussian', 'polynomial', 'contrast', 'filter', 'gradient', 'fft']
# METHODS = ['filter', 'fft']
METHODS = ['filter']

if CASE == 1:
    TITLE = 'All'
    TESTS = ['unagitated_1s',
             'unagitated_8s',
             'agitated_5volts_40mm_1s_new',
             'agitated_5volts_40mm_8s_new',
             'agitated_5volts_40mm_80s',
             'agitated_30volts_40mm_1s',
             'agitated_30volts_40mm_10s',
             'agitated_5volts_160mm_1s_new',
             'agitated_5volts_160mm_8s_new',
             'agitated_5volts_160mm_80s',
             'agitated_30volts_160mm_1s',
             'agitated_30volts_160mm_10s',
             'baseline']
    LABELS = ['unagitated 1s',
              'unagitated 8s',
              '0.1Hz 40mm 1s',
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

if CASE == 2:
    TITLE = 'Amptlitude vs Frequency'
    TESTS = ['unagitated_8s',
             'unagitated_1s',
             'agitated_5volts_40mm_8s_new',
             'agitated_5volts_160mm_8s_new',
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
    TITLE = 'New vs Old'
    TESTS = ['agitated_5volts_40mm_1s',
             'agitated_5volts_40mm_1s_new',
             'agitated_5volts_40mm_8s',
             'agitated_5volts_40mm_8s_new',
             'agitated_5volts_160mm_1s',
             'agitated_5volts_160mm_1s_new',
             'agitated_5volts_160mm_8s',
             'agitated_5volts_160mm_8s_new']
    LABELS = ['0.1Hz 40mm 1s',
              '0.1Hz 40mm 1s new',
              '0.1Hz 40mm 8s',
              '0.1Hz 40mm 8s new',
              '0.1Hz 160mm 1s',
              '0.1Hz 160mm 1s new',
              '0.1Hz 160mm 8s',
              '0.1Hz 160mm 8s new']

if CASE == 4:
    TITLE = 'Low Frequency'
    TESTS = ['unagitated_1s',
             'agitated_5volts_40mm_1s_new',
             'agitated_5volts_40mm_8s_new',
             'agitated_5volts_40mm_80s',
             'agitated_5volts_160mm_1s_new',
             'agitated_5volts_160mm_8s_new',
             'agitated_5volts_160mm_80s']
    LABELS = ['unagitated 1s',
              '0.1Hz 40mm 1s',
              '0.1Hz 40mm 8s',
              '0.1Hz 40mm 80s',
              '0.1Hz 160mm 1s',
              '0.1Hz 160mm 8s',
              '0.1Hz 160mm 80s']

if CASE == 5:
    TITLE = 'High Frequency'
    TESTS = ['unagitated_1s',
             'agitated_30volts_40mm_1s',
             'agitated_30volts_40mm_10s',
             'agitated_30volts_160mm_1s',
             'agitated_30volts_160mm_10s']
    LABELS = ['unagitated 1s',
              '1.0Hz 40mm 1s',
              '1.0Hz 40mm 10s',
              '1.0Hz 160mm 1s',
              '1.0Hz 160mm 10s']

if CASE == 6:
    TITLE = 'Low Amplitude'
    TESTS = ['unagitated_1s',
             'agitated_5volts_40mm_1s_new',
             'agitated_5volts_40mm_8s_new',
             'agitated_5volts_40mm_80s',
             'agitated_30volts_40mm_1s',
             'agitated_30volts_40mm_10s']
    LABELS = ['unagitated 1s',
              '0.1Hz 40mm 1s',
              '0.1Hz 40mm 8s',
              '0.1Hz 40mm 80s',
              '1.0Hz 40mm 1s',
              '1.0Hz 40mm 10s']

if CASE == 7:
    TITLE = 'High Amplitude'
    TESTS = ['unagitated_1s',
             'agitated_5volts_160mm_1s_new',
             'agitated_5volts_160mm_8s_new',
             'agitated_5volts_160mm_80s',
             'agitated_30volts_160mm_1s',
             'agitated_30volts_160mm_10s']
    LABELS = ['unagitated 1s',
              '0.1Hz 160mm 1s',
              '0.1Hz 160mm 8s',
              '0.1Hz 160mm 80s',
              '1.0Hz 160mm 1s',
              '1.0Hz 160mm 10s']

if CASE == 8:
    TITLE = 'Integration x10'
    TESTS = ['agitated_5volts_40mm_80s',
             'agitated_5volts_160mm_80s',
             'agitated_30volts_40mm_10s',
             'agitated_30volts_160mm_10s']
    LABELS = ['0.1Hz 40mm 80s',
              '0.1Hz 160mm 80s',
              '1.0Hz 40mm 10s',
              '1.0Hz 160mm 10s']

def main():
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

        for i, test in enumerate(TESTS):
            folder = FOLDER + test

            print cam, test
            new_object = NEW_OBJECTS or cam + '_obj.pkl' not in os.listdir(folder)
            if new_object:
                dark_folder = '../dark/'
                ambient_folder = '../ambient_1s/'
                if '1s_new' in folder:
                    ambient_folder = '../ambient_1s_new/'
                if '8s_new' in test:
                    ambient_folder = '../ambient_8s_new/'
                elif '8s' in test:
                    ambient_folder = '../ambient_8s/'
                if '10s' in test:
                    ambient_folder = '../ambient_10s/'
                if '80s' in test:
                    ambient_folder = '../ambient_80s/'
                save_new_object(folder, cam, ambient_folder, dark_folder)

            if NEW_DATA or new_object:
                set_new_data(folder, cam, methods,
                             fiber_method=FIBER_METHOD,
                             kernel_size=kernel)
        if 'fft' in methods:
            methods.remove('fft')
            save_fft_plot(FOLDER, TESTS, cam, LABELS, TITLE)

        save_modal_noise_bar_plot(FOLDER, TESTS, cam, LABELS, methods[0], TITLE)

if __name__ == '__main__':
    main()
