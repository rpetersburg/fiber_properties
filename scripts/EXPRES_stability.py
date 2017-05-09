import matplotlib.pyplot as plt
import numpy as np
import os
from fiber_properties import FiberImage, load_image_object, image_list
from functools import partial
from multiprocessing import Pool

plt.rc('font', size=32, family='serif')
plt.rc('figure', figsize=[20, 12.36])
plt.rc('xtick', labelsize=32)
plt.rc('ytick', labelsize=32)
plt.rc('text', usetex=True)
plt.rc('lines', lw=4)

PARALLELIZE = True
PROCESSES = 3

NUM_IMAGES = 150
CAMS = ['in', 'nf', 'ff']
# FOLDER = '../data/stability/2017-03-19 Stability Test/circular_200um/'
FOLDER = '../data/EXPRES/bare_octagonal/stability/'
# FOLDER = '../data/scrambling/2016-08-05 Prototype Core Extension 1/Shift_30/'
NF_METHOD = 'radius'
FF_METHOD = 'edge'
BIN_SIZE = 10

def plot_stability(data, cam):
    sigma = np.sqrt(np.std(data.x_diff)**2 + np.std(data.y_diff)**2)
    max_sg = np.median(data.diameter) / sigma

    plt.figure()

    plt.subplot(311)
    plt.plot(data.time, data.x_diff)
    plt.ylabel('x drift [um]')
    plt.title(r'$\sigma_{%s}= %.3f um, SG_{max} = %d$' % (cam, sigma, max_sg), fontsize=64)

    plt.subplot(312)
    plt.plot(data.time, data.y_diff)
    plt.ylabel('y drift [um]')

    plt.subplot(313)
    plt.plot(data.time, data.diameter)
    plt.xlabel('time [min]')
    plt.ylabel('diameter [um]')


def plot_stability_binned(data, cam, bin_size):
    plot_stability(data, cam)

    half_bin = bin_size / 2
    time = data.time[half_bin:len(data.time)-half_bin]
    x_diff = []
    y_diff = []
    diameter = []
    for i in xrange(half_bin, len(data.time)-half_bin):
        x_diff.append(np.array(data.x_diff[i-half_bin:i+half_bin]).mean())
        y_diff.append(np.array(data.y_diff[i-half_bin:i+half_bin]).mean())
        diameter.append(np.array(data.diameter[i-half_bin:i+half_bin]).mean())

    sigma = np.sqrt(np.std(x_diff)**2 + np.std(y_diff)**2)
    max_sg = np.median(diameter) / sigma
    print cam, 'max SG:', max_sg

    plt.subplot(311)
    plt.plot(time, x_diff, c='red')
    plt.title(r'$\sigma_{%s}= %.3f um, SG_{max} = %d$' % (cam, sigma, max_sg), fontsize=64)

    plt.subplot(312)
    plt.plot(time, y_diff, c='red')

    plt.subplot(313)
    plt.plot(time, diameter, c='red')

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
            obj = load_image_object(FOLDER + obj_file)
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
