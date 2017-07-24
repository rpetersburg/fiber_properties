from fiber_properties import (baseline_image, FiberImage,
                              plot_fft, save_plot, create_directory,
                              plot_modal_noise, join_paths, true_path)
import csv
import os
import numpy as np

DEFAULT_NUM = 10
DEFAULT_START = 0

def image_list(folder, cam, num=10, start=0, ext='.fit'):
    return [image_base(folder, cam, i) + ext for i in xrange(start, start+num, 1)]

def image_base(folder, cam, im):
    if folder and not folder.endswith('/'):
        folder += '/'
    return folder + cam + '_' + str(im).zfill(3)

def corrected_image_file(folder, cam, num=10, start=0, ext='.fit'):
    if folder and not folder.endswith('/'):
        folder += '/'
    front = image_base(folder, cam, start)
    back = '_corrected' + ext
    if start == 0 and not os.path.exists(image_base(folder, cam, num) + '.fit'):
        return folder + cam + back
    elif num == 1:
        return front + back
    else:
        return front + '-' + str(start+num-1).zfill(3) + back

def object_file(folder, cam, num=10, start=0):
    if folder and not folder.endswith('/'):
        folder += '/'
    front = image_base(folder, cam, start)
    back = '_obj.pkl'
    if start == 0 and not os.path.exists(image_base(folder, cam, num) + '.fit'):
        return folder + cam + back
    elif num == 1:
        return front + back
    else:
        return front + '-' + str(start+num-1).zfill(3) + back

def user_input(question='user input: ', valid_responses=['y', 'n']):
    response = raw_input(question)
    while response not in valid_responses:
        print 'invalid response'
        response = raw_input(question)
    return response

def save_new_object(folder, cam, ambient_folder=None, dark_folder=None,
                    overwrite='choose', num=10, start=0):
    if folder and not folder.endswith('/'):
        folder += '/'
    print 'saving object ' + object_file(folder, cam, num, start) + '...'
    images = image_list(folder, cam, num, start)

    ambient = None
    if ambient_folder:
        ambient = image_list(folder + ambient_folder, cam)
    dark = None
    if dark_folder:
        dark = image_list(folder + dark_folder, cam)

    response = 'y'
    if overwrite is not True and os.path.exists(object_file(folder, cam, num, start)):
        if overwrite is 'choose':
            response = user_input('overwrite old object? [y/n]: ')
        elif overwrite is False:
            response = 'n'
        else:
            raise RuntimeError('overwrite condition must be True or False')

    if response == 'y':
        im_obj = FiberImage(images, dark=dark, ambient=ambient, camera=cam)    
        im_obj.save_object(object_file(folder, cam, num, start))
        print 'object saved'
    else:
        print 'skipping'
    print

def set_new_data(folder, cam, methods, overwrite='choose', num=10, start=0, **kwargs):
    if folder and not folder.endswith('/'):
        folder += '/'
    im_obj = FiberImage(object_file(folder, cam, num, start))

    if 'rect' in folder and cam is 'nf':
        fiber_shape = 'rect'
    else:
        fiber_shape = 'circ'

    if isinstance(methods, basestring):
        methods = [methods]

    for method in methods:
        print 'setting ' + method + ' method...'

        response = 'y'
        if overwrite is not True and getattr(im_obj._modal_noise_info, method):
            if overwrite is 'choose':
                response = user_input('overwrite old data? [y/n]: ')
            elif overwrite is False:
                response = 'n'
            else:
                raise RuntimeError('overwrite condition must be True or False')

        if response == 'y':
            im_obj.set_modal_noise(method, fiber_shape=fiber_shape, **kwargs)
            im_obj.save_object(object_file(folder, cam, num, start))
            print method + ' method complete'
        else:
            print 'skipping'
    print

def save_new_image(folder, cam, num=10, start=0, ext='.fit'):
    print 'saving image ' + corrected_image_file(folder, cam, num, start, ext)
    im_obj = FiberImage(object_file(folder, cam, num, start))
    im_obj.save_image(corrected_image_file(folder, cam, ext))

def save_baseline_object(folder, cam, best_test, fiber_method='edge', kernel=None):
    print 'saving new baseline object'
    im_obj = FiberImage(object_file(folder + best_test, cam))
    baseline = baseline_image(im_obj, stdev=im_obj.get_dark_image().std(),
                              fiber_method=fiber_method, kernel_size=kernel)

    baseline_obj = FiberImage(baseline, camera=cam,
                              pixel_size=im_obj.pixel_size)
    save_new_image(folder, cam, ext='.fit')
    save_new_image(folder, cam, ext='.png')
    baseline_obj.save_object(object_file(folder, cam))

def save_fft_plot(folder, tests, cam, labels, title):
    print 'saving fft plot'
    fft_info_list = []
    for test in tests:
        im_obj = FiberImage(object_file(folder + test, cam))
        fft_info_list.append(im_obj.get_modal_noise(method='fft'))
    min_wavelength = im_obj.pixel_size / im_obj.magnification * 2.0
    max_wavelength = im_obj.get_fiber_radius(method='edge', units='microns')
    plot_fft(fft_info_list,
             labels=labels,
             min_wavelength=min_wavelength,
             max_wavelength=max_wavelength)
    save_plot(folder + 'analysis/' + title + ' ' + cam.upper() + ' FFT.png', dpi=600)
    # save_plot(folder + 'analysis/' + title + '/' + cam.upper() + '.pdf', dpi=600)

def save_modal_noise_data(folder, tests, cam, labels, methods, title=''):
    print 'saving modal noise data'
    modal_noise_info = [['cam', 'test'] + methods]

    for i, test in enumerate(tests):
        im_obj = FiberImage(object_file(folder + test, cam))
        modal_noise_info.append([cam, test])         
        for method in methods:
            modal_noise = im_obj.get_modal_noise(method)
            im_obj.save_object(object_file(folder + test, cam))
            print cam, test, method, modal_noise
            modal_noise_info[i+1].append(modal_noise)

    create_directory(folder + 'analysis/' + title + ' ' + cam.upper() + ' Data.csv')
    with open(folder + 'analysis/' + title + ' ' + cam.upper() + ' Data.csv', 'wb') as f:
        wr = csv.writer(f)
        wr.writerows(modal_noise_info)

def save_modal_noise_bar_plot(folder, tests, cam, bar_labels, method='filter', title='', labels=[''], num=1):
    modal_noise = []
    std = []
    for test in tests:
        mn = []
        for im in xrange(0, 10, num):
            im_obj = FiberImage(object_file(folder + test, cam, num, im))
            mn.append(im_obj.get_modal_noise(method=method))
        mn = np.array(mn)
        modal_noise.append(mn.mean())
        std.append(2.0*mn.std())
        # im_obj = FiberImage(object_file(folder + test, cam, 10, 0))
        # modal_noise.append(im_obj.get_modal_noise(method=method))
    plot_modal_noise([modal_noise], plot_type='bar', bar_labels=bar_labels, method=method, labels=labels, errors=[std])
    save_plot(folder + 'analysis/' + title + ' ' + cam.upper() + ' SNR.png')
    # save_plot(folder + 'analysis/' + title + ' ' + cam.upper() + ' SNR.pdf')

def save_modal_noise_line_plot(folder, tests, cam, labels=[''], method='filter', title=''):
    modal_noise = []
    for test in tests:
        mn = []
        for im in xrange(10):
            im_obj = FiberImage(object_file(folder + test, cam, im+1, 0))
            mn.append(im_obj.get_modal_noise(method=method))
        modal_noise.append(mn)
    plot_modal_noise(modal_noise, labels=labels, plot_type='line', method=method)
    save_plot(folder + 'analysis/' + title + ' ' + cam.upper() + ' SNR vs Time.png')
    # save_plot(folder + 'analysis/' + title + ' ' + cam.upper() + ' SNR vs Time.pdf')

def compress(data, selectors):
    return [d for d, s in zip(data, selectors) if s]

def save_modal_noise_inside(folder, cams, methods=['filter', 'fft'],
                            overwrite='choose', **kwargs):
    folder = true_path(folder) + '/'
    dir_list = os.listdir(folder)

    for item in dir_list:
        new_folder = join_paths(folder, item)
        if os.path.isdir(new_folder):
            save_modal_noise_inside(new_folder, cams, methods, overwrite, **kwargs)

    for cam in cams:
        data = compress(dir_list, [i.startswith(cam) and i.endswith('.fit')
                                   and 'corrected' not in i for i in dir_list])
        if data:
            max_num = max([int(i[-7:-4]) for i in data])
            if not any(cal_string in folder for cal_string in ['ambient', 'dark']):
                modal_noise = []
                fft_info_list = []
                modal_noise_time = []
                for i in xrange(max_num+1):
                    save_new_object(folder, cam, num=1, start=i, overwrite=overwrite)
                    set_new_data(folder, cam, methods, num=1, start=i,
                                 overwrite=overwrite, **kwargs)

                    save_new_object(folder, cam, num=i+1, start=0, overwrite=overwrite)
                    set_new_data(folder, cam, methods, num=i+1, start=0,
                                 overwrite=overwrite, **kwargs)

                    im_obj = FiberImage(object_file(folder, cam, num=1, start=i))
                    modal_noise.append(im_obj.get_modal_noise(method='filter'))
                    # fft_info_list.append(im_obj.get_modal_noise(method='fft'))

                    im_obj = FiberImage(object_file(folder, cam, num=i+1, start=0))
                    modal_noise_time.append(im_obj.get_modal_noise(method='filter'))

                labels = ['frame ' + str(i) for i in xrange(max_num+1)]
                plot_modal_noise([modal_noise], bar_labels=labels, plot_type='bar', method='filter')
                save_plot(folder + 'analysis/' + cam.upper() + ' SNR.png')
                plot_modal_noise([modal_noise_time], plot_type='line', method='filter')
                save_plot(folder + 'analysis/' + cam.upper() + ' SNR vs Time.png')

                # min_wavelength = im_obj.pixel_size / im_obj.magnification * 2.0
                # max_wavelength = im_obj.get_fiber_radius(method='edge', units='microns')
                # plot_fft(fft_info_list,
                #          labels=labels,
                #          min_wavelength=min_wavelength,
                #          max_wavelength=max_wavelength)
                # save_plot(folder + 'analysis/' + cam.upper() + ' FFT.png')

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Process script arguments')
    parser.add_argument('folder', help='folder location', type=str)
    parser.add_argument('methods', help='modal noise methods', nargs='*', default=['filter', 'fft'])
    parser.add_argument('-a', '--ambient', help='relative location of ambient folder', default=None)
    parser.add_argument('-d', '--dark', help='relative location of dark folder', default=None)
    parser.add_argument('-o', '--overwrite', help='overwrite all data', action='store_true')
    parser.add_argument('--new_object', help='create new objects', action='store_true')
    parser.add_argument('-i', '--image_ext', help='save corrected images with extension', default=None)
    parser.add_argument('-c', '--camera', help='cameras to use', nargs='*', type=str, default=['nf', 'ff'])
    parser.add_argument('-n', '--num', help='number of images per object', type=int, default=DEFAULT_NUM)
    parser.add_argument('-s', '--start', help='image number to start', type=int, default=DEFAULT_START)

    args = parser.parse_args()
    for cam in args.camera:
        print object_file('', cam, args.num, args.start)
        if args.new_object or (object_file('', cam, args.num, args.start)
                               not in os.listdir(args.folder)):
            save_new_object(folder=args.folder,
                            cam=cam,
                            ambient_folder=args.ambient,
                            dark_folder=args.dark,
                            start=args.start,
                            num=args.num)