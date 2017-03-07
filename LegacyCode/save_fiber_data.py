import re
import os
from fiber_properties import ImageAnalysis

def save_input_data(image, in_dark, in_ambient, in_flat): 
    file_name = 'Data/' + re.split('/|\.', image)[-2]

    if file_name + '_data.txt' in os.listdir('/'.join(image.split('/')[:-1])):
        print file_name + ' already saved'
        return

    obj = ImageAnalysis(image, in_dark, in_ambient, in_flat, camera='in')
    print file_name + ' initialized'

    obj.set_fiber_data(method='edge')
    obj.set_fiber_data(method='radius', tol=0.25, test_range=5)
    obj.set_fiber_data(method='gaussian')
    obj.set_fiber_centroid(method='full')

    obj.save_data(file_name=file_name)
    print file_name + ' saved, diameter:', obj.get_fiber_diameter(method='radius', units='microns')

def save_near_field_data(image, nf_dark, nf_ambient, nf_flat):
    file_name = 'Data/' + re.split('/|\.', image)[-2]

    if file_name + '_data.txt' in os.listdir('/'.join(image.split('/')[:-1])):
        print file_name + ' already saved'
        return

    obj = ImageAnalysis(image, nf_dark, nf_ambient, nf_flat camera='nf', threshold=500)
    print file_name + ' initialized'

    obj.set_fiber_data(method='edge')
    obj.set_fiber_data(method='radius', tol=0.25, test_range=5)
    obj.set_fiber_centroid(method='full')

    obj.save_data(file_name=file_name)
    print file_name + ' saved, diameter:', obj.get_fiber_diameter(method='radius', units='microns')

def save_far_field_data(image, ff_dark, ff_ambient):
    file_name = 'Data/' + re.split('/|\.', image)[-2]

    if file_name + '_data.txt' in os.listdir('/'.join(image.split('/')[:-1])):
        print file_name + ' already saved'
        return

    obj = ImageAnalysis(image, ff_dark, ff_ambient, camera='ff', magnification=1)
    print file_name + ' initialized'

    obj.set_fiber_data(method='gaussian')
    obj.set_fiber_centroid(method='full')

    obj.save_data(file_name=file_name)
    print file_name + ' saved, diameter:', obj.get_fiber_diameter(method='gaussian', units='microns')

if __name__ == '__main__':
    from input_output import image_list

    PARALLELIZE = True
    
    base_folder = 'Stability Measurements/2016-07-22 Stability Test/'
    ambient_folder = base_folder + 'Ambient/'
    dark_folder = base_folder + 'Dark/'
    flat_folder = base_folder + 'Flat/'
    folder = base_folder + 'stability_unagitated/'
    ext = '.fit'

    in_dark = image_list(dark_folder + 'in_')
    in_ambient = image_list(ambient_folder + 'in_')
    in_flat = image_list(flat_folder + 'in_', num=8)
    print 'IN calibration initialized'
    nf_dark = image_list(dark_folder + 'nf_')
    nf_ambient = image_list(ambient_folder + 'nf_')
    nf_flat = image_list(flat_folder + 'nf_', num=8)
    print 'NF calibration initialized'
    ff_dark = image_list(dark_folder + 'ff_')
    ff_ambient = image_list(ambient_folder + 'ff_0.001_')
    print 'FF calibration initialized'

    num_images = 100
    in_images = image_list(folder + 'in_', num=num_images)
    nf_images = image_list(folder + 'nf_', num=num_images)
    ff_images = image_list(folder + 'ff_', num=num_images)

    if PARALLELIZE:
        from multiprocessing import Pool
        from functools import partial
        pool = Pool(processes=6)
        pool.map(partial(save_input_data, in_dark=in_dark, in_ambient=in_ambient, in_flat=in_flat), in_images)
        print 'IN data initialized'
        pool.map(partial(save_near_field_data, nf_dark=nf_dark, nf_ambient=nf_ambient, nf_flat=nf_flat), nf_images)
        print 'NF data initialized'
        pool.map(partial(save_far_field_data, ff_dark=ff_dark, ff_ambient=ff_ambient), ff_images)
        print 'FF data initialized'

    else:    
        for image in in_images:
            save_input_data(image, in_dark, in_ambient, in_flat)
        print 'IN data initialized'

        for image in nf_images:
            save_near_field_data(image, nf_dark, nf_ambient, nf_flat)
        print 'NF data initialized'

        for image in ff_images:
            save_far_field_data(image, ff_dark, ff_ambient)
        print 'FF data initialized'