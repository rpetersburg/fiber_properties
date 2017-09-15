from modal_noise_script import save_modal_noise_inside
import os

FOLDER = '../data/modal_noise/rec_fiber_amp_tests/'
CAMERAS = ['nf']
METHODS = ['filter']
FIBER_METHOD = 'full'
KERNEL = 51
OVERWRITE = False

def main():
    folder = FOLDER
    cams = CAMERAS
    methods = METHODS
    kernel = KERNEL
    fiber_method = FIBER_METHOD
    save_modal_noise_inside(folder, cams, methods,
                            overwrite=OVERWRITE,
                            kernel_size=kernel,
                            fiber_method=fiber_method,
                            ambient_folder='../ambient_2s/',
                            dark_folder=None)

if __name__ == '__main__':
    main()