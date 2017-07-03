from modal_noise_script import save_new_object, set_new_data
from fiber_properties import join_paths, true_path
import os

FOLDER = '../data/modal_noise/Kris_data/circular_100um/'
CAMERAS = ['nf', 'ff']
METHODS = ['filter', 'fft']
FIBER_METHOD = 'edge'
KERNEL = 51
OVERWRITE = False

def compress(data, selectors):
    return [d for d, s in zip(data, selectors) if s]

def save_modal_noise_inside(folder, cams, methods, **kwargs):
    folder = true_path(folder)
    dir_list = os.listdir(folder)

    for item in dir_list:
        new_folder = join_paths(folder, item)
        if os.path.isdir(new_folder):
            save_modal_noise_inside(new_folder, cams, methods, **kwargs)

    for cam in cams:
        data = compress(dir_list, [i.startswith(cam) and i.endswith('.fit')
                                   and 'corrected' not in i for i in dir_list])
        if data:
            max_num = max([int(i[-7:-4]) for i in data])
            if any(cal_string in folder for cal_string in ['ambient', 'dark']):
                for i in xrange(max_num+1):
                    save_new_object(folder, cam, num=1, start=i, overwrite=OVERWRITE)
            else:
                for i in xrange(max_num+1):
                    save_new_object(folder, cam, num=1, start=i, overwrite=OVERWRITE)
                    set_new_data(folder, cam, methods, num=1, start=i, overwrite=OVERWRITE, **kwargs)
                save_new_object(folder, cam, num=max_num+1, start=0, overwrite=OVERWRITE)
                set_new_data(folder, cam, methods, num=max_num+1, start=0, overwrite=OVERWRITE, **kwargs)

def main():
    folder = FOLDER
    cams = CAMERAS
    methods = METHODS
    kernel = KERNEL
    fiber_method = FIBER_METHOD

    save_modal_noise_inside(folder, cams, methods,
                            kernel_size=kernel,
                            fiber_method=fiber_method)

if __name__ == '__main__':
    main()