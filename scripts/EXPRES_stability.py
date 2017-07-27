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

NEW_DATA = False
PARALLELIZE = True
PROCESSES = 3

NUM_IMAGES = 150
CAMS = ['nf']
# FOLDER = '../data/stability/2017-03-19 Stability Test/circular_200um/'
FOLDER = '../data/EXPRES/bare_octagonal/stability/'
# FOLDER = '../data/scrambling/2016-08-05 Prototype Core Extension 1/Shift_30/'
NF_METHOD = 'full'
FF_METHOD = 'full'
BIN_SIZE = 10
OPTIONS = {'units': 'microns',
           'threshold': 1000,
           'kernel': 9,
           'radius_tol': .03,
           'radius_range': 64,
           'center_tol': .03,
           'center_range': 64}

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
        obj = FiberImage(im_file)
        obj.set_fiber_center(method=method, **OPTIONS)
        obj.set_fiber_centroid(method=method, **OPTIONS)
        obj.save_object(FOLDER + obj_file)

if __name__ == "__main__":
    data = {}

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
            data[cam].center.append(obj.get_fiber_center(method=method, **OPTIONS))
            data[cam].centroid.append(obj.get_fiber_centroid(method=method, **OPTIONS))
            data[cam].x_diff.append(data[cam].centroid[-1].x - data[cam].center[-1].x)
            data[cam].y_diff.append(data[cam].centroid[-1].y - data[cam].center[-1].y)
            # data[cam].x_diff.append(data[cam].centroid[-1].x)
            # data[cam].y_diff.append(data[cam].centroid[-1].y)
            data[cam].diameter.append(obj.get_fiber_diameter(method=method, **OPTIONS))
            data[cam].time.append(obj.date_time)
            obj.save_object(FOLDER + obj_file)

        avg_x_diff = np.median(data[cam].x_diff)
        avg_y_diff = np.median(data[cam].y_diff)
        init_time = np.copy(data[cam].time[0])
        for i in xrange(NUM_IMAGES):
            data[cam].x_diff[i] -= avg_x_diff
            data[cam].y_diff[i] -= avg_y_diff
            data[cam].time[i] -= init_time
            data[cam].time[i] = data[cam].time[i].total_seconds() / 60.0
        plot_stability(data[cam], cam)
        plt.savefig(FOLDER + cam + '_stability.png')
        plot_stability_binned(data[cam], cam, BIN_SIZE)
        plt.savefig(FOLDER + cam + '_stability_binned.png')
