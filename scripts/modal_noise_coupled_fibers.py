from modal_noise_script import (save_new_object, set_new_data,
                                save_baseline_object, save_fft_plot,
                                save_modal_noise_data)
from copy import deepcopy
import os

NEW_DATA = True
NEW_OBJECTS = False
NEW_BASELINE = False
FIBER_METHOD = 'edge'
CAMERAS = ['ff']
CASE = 3
KERNEL = 51
METHODS = ['tophat', 'gaussian', 'polynomial', 'contrast', 'filter', 'gradient', 'fft']
FOLDER = "C:/Libraries/Box Sync/ExoLab/Fiber_Characterization/Image Analysis/data/modal_noise/coupled_fibers/"

if CASE == 1:
    TITLE = 'Modal Noise 200-200um'
    FOLDER += '200-200um_test2/'
    TESTS = ['unagitated',
             'agitated_first',
             'agitated_second',
             'agitated_both',
             'baseline']
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
    LABELS = ['unagitated',
              'circular agitated',
              'octagonal agitated',
              'both agitated',
              'baseline']

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

      base_i = None
      for i, test in enumerate(TESTS):
         if 'baseline' in test:
             base_i = i
             continue       

         print cam, test
         new_object = NEW_OBJECTS or cam + '_obj.pkl' not in os.listdir(FOLDER + test + '/')
         if new_object:
            ambient = 'ambient/'
            dark = 'dark/'
            save_new_object(FOLDER + test, cam, ambient, dark)

         if NEW_DATA or new_object:
             set_new_data(FOLDER + test, cam, methods,
                          fiber_method=FIBER_METHOD,
                          kernel_size=kernel)

      if 'fft' in methods:
         methods.remove('fft')
         save_fft_plot(FOLDER, TESTS, cam, LABELS, TITLE)

      save_modal_noise_data(FOLDER, TESTS, cam, LABELS, methods, TITLE)


if __name__ == '__main__':
   main()
