from fiber_properties import (image_list, baseline_image,
                              FiberImage, plot_fft, save_plot)
import csv

def image_file(folder, test, cam):
    return folder + test + '/' + cam + '_corrected.fit'

def object_file(folder, test, cam):
    return folder + test + '/' + cam + '_obj.pkl'

def save_new_object(folder, test, cam, ambient_folder='ambient/', dark_folder='dark/'):
    print 'saving new object'
    images = image_list(folder + test + '/' + cam + '_')

    ambient = image_list(folder + ambient_folder + cam + '_')
    dark = image_list(folder + dark_folder + cam + '_')

    im_obj = FiberImage(images, dark=dark, ambient=ambient, camera=cam)
    im_obj.save_object(object_file(folder, test, cam))

def set_new_data(folder, test, cam, methods, fiber_method='edge', kernel=None):
    if cam == 'ff':
        kernel = None
    im_obj = FiberImage(object_file(folder, test, cam))
    radius_factor = None
    if 'rectang' in test:
        radius_factor = 0.3
        
    print 'setting new data'
    for method in methods:
        print 'setting method ' + method
        im_obj.set_modal_noise(method, fiber_method=fiber_method,
                               kernel_size=kernel, radius_factor=radius_factor)
    im_obj.save_object(object_file(folder, test, cam))
    im_obj.save_image(image_file(folder, test, cam))
    im_obj.save_image(image_file(folder, test, cam)[:-3] + 'png')
    print

def save_baseline_object(folder, test, cam, best_test, fiber_method='edge', kernel=None):
    print 'saving new baseline object'
    im_obj = FiberImage(object_file(folder, best_test, cam))
    baseline = baseline_image(im_obj, stdev=im_obj.get_dark_image().std(),
                              fiber_method=fiber_method, kernel_size=kernel)

    baseline_obj = FiberImage(baseline, camera=cam,
                              pixel_size=im_obj.pixel_size)
    baseline_obj.save_image(image_file(folder, test, cam))
    baseline_obj.save_image(image_file(folder, test, cam)[:-3] + 'png')
    baseline_obj.save_object(object_file(folder, test, cam))

def save_fft_plot(folder, tests, cam, labels, title):
    print 'saving fft plot'
    fft_info_list = []
    for test in tests:
        im_obj = FiberImage(object_file(folder, test, cam))
        fft_info_list.append(im_obj.get_modal_noise(method='fft'))
    min_wavelength = im_obj.pixel_size / im_obj.magnification * 2.0
    max_wavelength = im_obj.get_fiber_radius(method='edge', units='microns')
    plot_fft(fft_info_list,
             labels=labels,
             min_wavelength=min_wavelength,
             max_wavelength=max_wavelength)
    save_plot(folder + 'analysis/' + title + '/' + cam.upper() + '.png', dpi=600)
    save_plot(folder + 'analysis/' + title + '/' + cam.upper() + '.pdf', dpi=600)

def save_modal_noise_data(folder, tests, cam, methods, title):
    print 'saving modal noise data'
    modal_noise_info = [['cam', 'test'] + methods]
    for i, test in enumerate(tests):
        im_obj = FiberImage(object_file(folder, test, cam))
        modal_noise_info.append([cam, test])         
        for method in methods:
            modal_noise = im_obj.get_modal_noise(method)
            im_obj.save_object(object_file(folder, test, cam))
            print cam, test, method, modal_noise
            modal_noise_info[i+1].append(modal_noise)

    with open(folder + 'analysis/' + title + '/' + cam.upper() + ' Data.csv', 'wb') as f:
        wr = csv.writer(f)
        wr.writerows(modal_noise_info)