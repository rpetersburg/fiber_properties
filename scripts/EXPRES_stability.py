import matplotlib.pyplot as plt
import numpy as np
import os
from fiber_properties import ImageAnalysis, load_image_object, image_list

plt.rc('font', size=16, family='serif')
plt.rc('figure', figsize=[20, 12.36])
plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('lines', lw=4)

NUM_IMAGES = 500
CAMS = ['nf']
FOLDER = '../data/stability/2017-03-19 Stability Test/circular_200um/'
NF_METHOD = 'radius'
FF_METHOD = 'gaussian'
BIN_SIZE = 10

def plot_stability(data, cam):
    plt.figure() 
    plt.title(cam + 'stability')
    plt.subplot(311)
    plt.plot(data.time, data.x_diff)
    plt.ylabel('X Centroid Shift [um]')
    plt.title('std: ' + str(np.std(data.x_diff)))

    plt.subplot(312)
    plt.plot(data.time, data.y_diff)
    plt.ylabel('Y Centroid Shift [um]')
    plt.title('std: ' + str(np.std(data.y_diff)))

    plt.subplot(313)
    plt.plot(data.time, data.diameter)
    plt.xlabel('Time [min]')
    plt.ylabel('Diameter Shift [um]')
    plt.title('std: ' + str(np.std(data.diameter)))

    plt.suptitle(cam + ' stability')

def plot_stability_binned(data, cam, bin_size):
    plot_stability(data, cam)

    time = data.time[:len(data.time)-bin_size]
    x_diff = []
    y_diff = []
    diameter = []
    for i in xrange(0, len(data.time)-bin_size):
        x_diff.append(np.array(data.x_diff[i:i+bin_size]).mean())
        y_diff.append(np.array(data.y_diff[i:i+bin_size]).mean())
        diameter.append(np.array(data.diameter[i:i+bin_size]).mean())

    plt.subplot(311)
    plt.plot(time, x_diff, c='red')
    plt.title('std: ' + str(np.array(x_diff).std()))

    plt.subplot(312)
    plt.plot(time, y_diff, c='red')
    plt.title('std: ' + str(np.array(y_diff).std()))

    plt.subplot(313)
    plt.plot(time, diameter, c='red')
    plt.title('std:  ' + str(np.array(diameter).std()))

class StabilityInfo(object):
    def __init__(self):
        self.centroid = []
        self.center = []
        self.x_diff = []
        self.y_diff  = []
        self.diameter = []
        self.time = []

if __name__ == "__main__":
    data = {}
    data['spot'] = StabilityInfo()

    for cam in CAMS:
        data[cam] = StabilityInfo()
        if cam == 'in' or cam == 'nf':
            method = NF_METHOD
        else:
            method = FF_METHOD

        for i in xrange(NUM_IMAGES):
            obj_file = cam + '_obj_' + str(i).zfill(3) + '.pkl'
            if obj_file not in os.listdir(FOLDER):
                print 'saving ' + cam + '_' + str(i)
                im_file = FOLDER + cam + '_' + str(i).zfill(3) + '.fit'
                ImageAnalysis(im_file, threshold=1000).save_object(FOLDER + obj_file)

            print 'loading ' + cam + '_' + str(i)
            obj = load_image_object(FOLDER + obj_file)
            data[cam].center.append(np.array(obj.get_fiber_center(method=method, units='microns')))
            data[cam].centroid.append(np.array(obj.get_fiber_centroid(method=method, units='microns')))
            data[cam].x_diff.append(data[cam].centroid[-1][1] - data[cam].center[-1][1])
            data[cam].y_diff.append(data[cam].centroid[-1][0] - data[cam].center[-1][0])
            # data[cam].x_diff.append(data[cam].center[-1][1])
            # data[cam].y_diff.append(data[cam].center[-1][0])
            data[cam].diameter.append(obj.get_fiber_diameter(method=method, units='microns'))
            data[cam].time.append(obj.get_image_info('date_time'))
            obj.save_object()
            if cam == 'in':
                data['spot'].center.append(np.array(obj.get_fiber_center(method='gaussian', units='microns')))
                data['spot'].centroid.append(np.array(obj.get_fiber_centroid(method='gaussian', units='microns')))
                data['spot'].x_diff.append(data['spot'].center[-1][1] - data[cam].center[-1][1])
                data['spot'].y_diff.append(data['spot'].center[-1][0] - data[cam].center[-1][0])
                # data['spot'].x_diff.append(data['spot'].center[i][1])
                # data['spot'].y_diff.append(data['spot'].center[i][0])
                data['spot'].diameter.append(obj.get_fiber_diameter(method='gaussian', units='microns'))
                data['spot'].time.append(obj.get_image_info('date_time'))
                obj.save_object()

    if 'in' in CAMS:
        CAMS += ['spot']

    for cam in CAMS:
        init_center = np.copy(data[cam].center[0])
        init_centroid = np.copy(data[cam].centroid[0])
        init_x_diff = np.copy(data[cam].x_diff[0])
        init_y_diff = np.copy(data[cam].y_diff[0])
        init_diameter = np.copy(data[cam].diameter[0])
        init_time = np.copy(data[cam].time[0])
        for i in xrange(NUM_IMAGES):
            data[cam].center[i] -= init_center
            data[cam].centroid[i] -= init_centroid
            data[cam].x_diff[i] -= init_x_diff
            data[cam].y_diff[i] -= init_y_diff
            data[cam].diameter[i] -= init_diameter
            data[cam].time[i] -= init_time
            data[cam].time[i] = data[cam].time[i].total_seconds() / 60.0
        plot_stability(data[cam], cam)
        plt.savefig(FOLDER + cam + '_stability.png')
        plot_stability_binned(data[cam], cam, BIN_SIZE)
        plt.savefig(FOLDER + cam + '_stability_binned.png')
