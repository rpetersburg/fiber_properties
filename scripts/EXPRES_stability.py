import matplotlib.pyplot as plt
import numpy as np
import os
from fiber_properties import FiberImage, image_list, plot_stability, plot_stability_binned
from functools import partial
from multiprocessing import Pool

plt.rc('font', size=32, family='serif')
plt.rc('figure', figsize=[18,18])
plt.rc('xtick', labelsize=32)
plt.rc('ytick', labelsize=32)
plt.rc('text', usetex=True)
plt.rc('lines', lw=4)

NEW_DATA = True
PARALLELIZE = True
PROCESSES = 3

NUM_IMAGES = 150
CAMS = ['nf']
# FOLDER = '../data/stability/2017-03-19 Stability Test/circular_200um/'
FOLDER = '../data/EXPRES/bare_octagonal/stability/'
# FOLDER = '../data/scrambling/2016-08-05 Prototype Core Extension 1/Shift_30/'
NF_METHOD = 'full'
FF_METHOD = 'edge'
BIN_SIZE = 10

class StabilityInfo(object):
    def __init__(self):
        self.centroid = []
        self.center = []
        self.x_diff = []
        self.y_diff  = []
        self.diameter = []
        self.time = []

def save_objects(i, cam, method):
    obj_file = cam + '_obj_' + str(i).zfill(3) + '.pkl'
    if obj_file not in os.listdir(FOLDER):
        print 'saving ' + cam + '_' + str(i).zfill(3)
        im_file = FOLDER + cam + '_' + str(i).zfill(3) + '.fit'
        obj = FiberImage(im_file, threshold=1000)
        obj.set_fiber_center(method=method, 
                             radius_tol=.03, radius_range=64,
                             center_tol=.03, center_range=64)
        obj.set_fiber_centroid(method=method)
        obj.save_object(FOLDER + obj_file)

if __name__ == "__main__":
    data = {}
    data['spot'] = StabilityInfo()

    for cam in CAMS:
        data[cam] = StabilityInfo()
        if cam == 'in' or cam == 'nf':
            method = NF_METHOD
        else:
            method = FF_METHOD

        if NEW_DATA:
            if PARALLELIZE:
                pool = Pool(processes=PROCESSES)
                pool.map(partial(save_objects, cam=cam, method=method),
                         range(NUM_IMAGES))
            else:
                for i in xrange(NUM_IMAGES):
                    save_objects(i, cam, method)

        for i in xrange(NUM_IMAGES):
            obj_file = cam + '_obj_' + str(i).zfill(3) + '.pkl'

            print 'loading ' + cam + '_' + str(i).zfill(3)
            obj = FiberImage(FOLDER + obj_file)
            data[cam].center.append(obj.get_fiber_center(method=method, units='microns'))
            data[cam].centroid.append(obj.get_fiber_centroid(method=method, units='microns'))
            data[cam].x_diff.append(data[cam].centroid[-1].x - data[cam].center[-1].x)
            data[cam].y_diff.append(data[cam].centroid[-1].y - data[cam].center[-1].y)
            data[cam].diameter.append(obj.get_fiber_diameter(method=method, units='microns'))
            data[cam].time.append(obj.date_time)
            obj.save_object(FOLDER + obj_file)
            if cam == 'in':
                data['spot'].center.append(obj.get_fiber_center(method='gaussian', units='microns'))
                data['spot'].centroid.append(obj.get_fiber_centroid(method='gaussian', units='microns'))
                data['spot'].x_diff.append(data['spot'].center[-1].x - data[cam].center[-1].x)
                data['spot'].y_diff.append(data['spot'].center[-1].y - data[cam].center[-1].y)
                data['spot'].diameter.append(obj.get_fiber_diameter(method='gaussian', units='microns'))
                data['spot'].time.append(obj.date_time)
                obj.save_object(FOLDER + obj_file)

    if 'in' in CAMS:
        CAMS += ['spot']

    for cam in CAMS:
        init_x_diff = np.copy(data[cam].x_diff[0])
        init_y_diff = np.copy(data[cam].y_diff[0])
        init_time = np.copy(data[cam].time[0])
        for i in xrange(NUM_IMAGES):
            data[cam].x_diff[i] -= init_x_diff
            data[cam].y_diff[i] -= init_y_diff
            data[cam].time[i] -= init_time
            data[cam].time[i] = data[cam].time[i].total_seconds() / 60.0
        plot_stability(data[cam], cam)
        plt.savefig(FOLDER + cam + '_stability.png')
        plot_stability_binned(data[cam], cam, BIN_SIZE)
        plt.savefig(FOLDER + cam + '_stability_binned.png')
