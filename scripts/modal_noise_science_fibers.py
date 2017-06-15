from modal_noise_script import (save_new_object, set_new_data,
                                save_baseline_object, save_fft_plot,
                                save_modal_noise_data)
from copy import deepcopy
import os

NEW_DATA = True
NEW_OBJECTS = True
NEW_BASELINE = True
FIBER_METHOD = 'edge'
KERNEL = 31
CAMERAS = ['nf']
CASE = 1
METHODS = ['tophat', 'gaussian', 'polynomial', 'contrast', 'filter', 'gradient', 'fft']
FOLDER = "C:/Libraries/Box Sync/ExoLab/Fiber_Characterization/Image Analysis/data/modal_noise/geometries/"

if CASE == 1:
    TITLE = 'Modal Noise Science Octagonal'
    TESTS = ['oct_60_unagitated',
             'oct_60_agitated',
             'oct_baseline']
if CASE == 2:
    TITLE = 'Modal Noise Science Rectangular'
    TESTS = ['rect_33x132_unagitated',
             'rect_33x132_agitated',
             'rect_33x132_baseline']
if CASE in [1,2]:
    LABELS = ['unagitated',
              'agitated',
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
            if 'baseline' in test:
                base_i = i
                continue       

            print cam, test
            new_object = NEW_OBJECTS or cam + '_obj.pkl' not in os.listdir(FOLDER + test + '/')
            if new_object:
                save_new_object(FOLDER, test, cam)

            if NEW_DATA or new_object:
                set_new_data(FOLDER, test, cam, methods, FIBER_METHOD, KERNEL)

        if base_i is not None:
            new_baseline = NEW_BASELINE or cam + '_obj.pkl' not in os.listdir(FOLDER + TESTS[base_i] + '/')
            if new_baseline:
                save_baseline_object(FOLDER, TESTS[base_i], cam, TESTS[base_i-1], FIBER_METHOD, KERNEL)

            if NEW_DATA or new_baseline:
                set_new_data(FOLDER, TESTS[base_i], cam, methods, FIBER_METHOD, KERNEL)

        if 'fft' in methods:
            methods.remove('fft')
            save_fft_plot(FOLDER, TESTS, cam, LABELS, TITLE)

        save_modal_noise_data(FOLDER, TESTS, cam, methods, TITLE)
