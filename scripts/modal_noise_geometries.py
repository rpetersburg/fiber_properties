from modal_noise_script import (save_new_object, set_new_data,
                                save_baseline_object, save_fft_plot,
                                save_modal_noise_data)
from copy import deepcopy
import os

NEW_DATA = True
NEW_OBJECTS = False
NEW_BASELINE = False
FIBER_METHOD = 'edge'
CAMERAS = ['nf']
CASE = 3
METHODS = ['tophat', 'gaussian', 'polynomial', 'contrast', 'filter', 'gradient', 'fft']
FOLDER = "C:/Libraries/Box Sync/ExoLab/Fiber_Characterization/Image Analysis/data/modal_noise/Kris_data/"

if CASE == 1:
    TITLE = 'Modal Noise Octagonal 100um'
    FOLDER += "octagonal_100um/"
if CASE == 2:
    TITLE = 'Modal Noise Circular 100um'
    FOLDER += "circular_100um/"
if CASE == 3:
    TITLE = 'Modal Noise Rectangular 100x300um'
    FOLDER += "rectangular_100x300um/"
if CASE == 7:
    TITLE = 'Modal Noise Octagonal 200um'
    FOLDER += "octagonal_200um/"
if CASE == 8:
    TITLE = 'Modal Noise Circular 200um'
    FOLDER += "circular_200um/"
if CASE in [1,2,3,7,8]:
    TESTS = ['unagitated',
             'linear_agitation',
             'circular_agitation',
             'coupled_agitation',
             'baseline']
    LABELS = ['unagitated',
              'linear agitation',
              'circular agitation',
              'coupled agitation',
              'baseline']
if CASE == 4:
    TITLE = 'Modal Noise Unagitated'
    TESTS = ['octagonal_200um/unagitated',
             'circular_200um/unagitated',
             'rectangular_100x300um/unagitated']
if CASE == 5:
    TITLE = 'Modal Noise Coupled Agitation'
    TESTS = ['octagonal_200um/coupled_agitation',
             'circular_200um/coupled_agitation',
             'rectangular_100x300um/coupled_agitation']
if CASE == 6:
    TITLE = 'Modal Noise Geometry Baselines'
    TESTS = ['octagonal_200um/baseline',
             'circular_200um/baseline',
             'rectangular_100x300um/baseline']
if CASE in [4,5,6]:
    LABELS = ['octagonal',
              'circular',
              'rectangular']

if CASE in [4,5]:
    CAL = ['octagonal_200um/',
           'circular_200um/',
           'rectangular_100x300um/']
else:
    CAL = ['' for i in xrange(len(TESTS))]

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
                save_new_object(FOLDER, test, cam, ambient_folder=CAL[i]+'ambient/', dark_folder=CAL[i]+'dark/')

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
            save_fft_plot(FOLDER, TESTS, cam, LABELS, TITLE)

        save_modal_noise_data(FOLDER, TESTS, cam, methods, TITLE)
