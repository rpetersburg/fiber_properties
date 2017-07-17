from modal_noise_script import save_modal_noise_inside
import os

FOLDER = '../data/modal_noise/Kris_data/'
CAMERAS = ['nf', 'ff']
METHODS = ['filter', 'fft']
FIBER_METHOD = 'edge'
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
                            fiber_method=fiber_method)

if __name__ == '__main__':
    main()