from modal_noise_script import (save_new_object, set_new_data,
                                save_baseline_object, save_fft_plot,
                                save_modal_noise_data)
from copy import deepcopy
import os

NEW_DATA = False
NEW_OBJECTS = False
NEW_BASELINE = False
FIBER_METHOD = 'edge'
CAMERAS = ['nf', 'ff']
CASE = 3
KERNEL = 51
METHODS = ['filter', 'fft']
FOLDER = "../data/modal_noise/coupled_fibers/"

if CASE == 1:
    TITLE = 'Modal Noise 200-200um'
    FOLDER += '200-200um_test2/'
    TESTS = ['unagitated',
             'agitated_first',
             'agitated_second',
             'agitated_both',
             'baseline']
    AMBIENTS = ['../ambient/' for i in len(TESTS)]
    DARKS = ['../dark/' for i in len(TESTS)]
    LABELS = ['unagitated',
              'first agitated',
              'second agitated',
              'both agitated',
              'baseline']
if CASE == 2:
    TITLE = 'Modal Noise 100-200um'
    FOLDER += '100-200um/'
    TESTS = ['unagitated',
             'agitated_first_100um',
             'agitated_second_200um',
             'agitated_both',
             'baseline']
    AMBIENTS = ['../ambient/' for i in len(TESTS)]
    DARKS = ['../dark/' for i in len(TESTS)]
    LABELS = ['unagitated',
              '100um agitated',
              '200um agitated',
              'both agitated',
              'baseline']
if CASE == 3:
    TITLE = 'Modal Noise Oct-Circ 200um'
    FOLDER += 'oct-circ-200um/'
    TESTS = ['unagitated',
             'agitated_circ',
             'agitated_oct',
             'agitated_both',
             'baseline']
    AMBIENTS = ['../ambient/' for i in len(TESTS)]
    DARKS = ['../dark/' for i in len(TESTS)]
    LABELS = ['unagitated',
              'circular agitated',
              'octagonal agitated',
              'both agitated',
              'baseline']

if CASE == 4:
    TITLE = 'Coupled Fiber'
    FOLDER += 'coupled_fibers/'
    TITLES = ['200um-200um',
              '100um-200um',
              'Oct-Circ',
              '200um']
    FOLDERS = [FOLDER + '200-200um_test2/',
               FOLDER + '100-200um/',
               FOLDER + 'oct-circ-200um/',
               FOLDER + '../Kris_data/circular_200um/']
    TESTS = [['agitated_first', 'agitated_second', 'agitated_both'],
             ['agitated_first_100um', 'agitated_second_200um', 'agitated_both'],
             ['agitated_oct', 'agitated_circ', 'agitated_both'],
             ['linear_agitation' for i in xrange(3)]]
    LABELS = ['first agitated', 'second agitated', 'both agitated']

def main():
    print TITLE
    print
    for cam in CAMERAS:

        kernel = KERNEL

        for test, ambient, dark in zip(TESTS, AMBIENTS, DARKS):
            folder = FOLDER + test
            print cam, test
            save_modal_noise_inside(folder, [cam], METHODS, overwrite=OVERWRITE,
                                    ambient_folder=ambient, dark_folder=dark,
                                    kernel_size=KERNEL, fiber_method=FIBER_METHOD)
            # new_object = NEW_OBJECTS or cam + '_obj.pkl' not in os.listdir(FOLDER + test + '/')
            # if new_object:
            #     ambient = 'ambient/'
            #     dark = 'dark/'
            #     save_new_object(FOLDER + test, cam, ambient, dark)

            # if NEW_DATA or new_object:
            #      set_new_data(FOLDER + test, cam, methods,
            #                   fiber_method=FIBER_METHOD,
            #                   kernel_size=kernel)

        if 'fft' in METHODS:
            save_fft_plot(FOLDER, TESTS, cam, LABELS, TITLE)
        else:
            save_modal_noise_bar_plot(FOLDER, TESTS, cam, LABELS, , TITLE + ' 1x', num=1)
            save_modal_noise_bar_plot(FOLDER, TESTS, cam, LABELS, method, TITLE + ' 10x', num=10)
            save_modal_noise_line_plot(FOLDER, TESTS, cam, LABELS, method, TITLE)


if __name__ == '__main__':
   main()
