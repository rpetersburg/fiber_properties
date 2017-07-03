from fiber_properties import FiberImage
import matplotlib.pyplot as plt
import os
import numpy as np
from multiprocessing import Pool

PLOT_FIBER_CENTROID = True
NEW_DATA = False
PROCESSES = 4
NUM_IMAGES = 1
FOLDER = '/home/deggerman/data/modal_noise/rv_error/'
CAMERAS = ['nf']
METHOD = 'circle'
CENTER_RANGE = 10  # for circle method. Default None

if CASE == 1:
    FOLDER += 'coupled_agitation/'
if CASE == 2:
    FOLDER += 'LED/'
if CASE == 3:
    FOLDER += 'slow_agitation/'

def multi(processes=PROCESSES, camera=CAMERAS, num_images=NUM_IMAGES):
    p = Pool(processes)
    for cam in camera:
        objects = []
        for i in xrange(0, 100, num_images):
            object_file = cam + '_' + str(i).zfill(3) + '_obj.pkl'
            objects.append(object_file)
        p.map(find_fiber_center, objects)

def find_fiber_center(object_file, folder=FOLDER, cameras=CAMERAS, num_images=NUM_IMAGES, new_data=NEW_DATA, meth=METHOD, plot_fiber_centroid=PLOT_FIBER_CENTROID, cent_range=CENTER_RANGE):

    cam = object_file[:2]
    ject = object_file[:6]
    if object_file not in os.listdir(folder) or new_data:
        images = [folder + ject + '.fit']
        ambient = [folder + 'ambient/' + cam + '_' + str(j).zfill(3) + '.fit' for j in xrange(10)]
        dark = [folder + 'dark/' + cam + '_' + str(j).zfill(3) + '.fit' for j in xrange(10)]
        im_obj = FiberImage(images, ambient=ambient, dark=dark, camera=cam)
        im_obj.save_object(folder + object_file)

    object_file = folder + object_file
    im_obj = FiberImage(object_file)

    print('Getting fiber center for image %s...' % object_file[-14:])
    fiber_centroid = im_obj.get_fiber_center(method=meth, units='microns', center_range=cent_range) - im_obj.get_fiber_centroid(method=meth, units='microns', center_range=cent_range)
    im_obj.save_object(object_file)


if __name__ == '__main__':
    multi()
