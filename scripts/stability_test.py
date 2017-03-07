import matplotlib.pyplot as plt
from copy import deepcopy
import numpy as np
from fiber_properties import ImageAnalysis, convert_pixels_to_microns
from datetime import datetime
import os

plt.rc('font', size=16, family='serif')
plt.rc('figure', figsize=[20, 12.36])
plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('lines', lw=4)

NUM_IMAGES = 100
SAVE_FOLDER = 'Stability Measurements/Bench '
TESTS = ['bench agitated', 'bench unagitated']
CAMERAS = ['nf', 'ff']

FOLDER = {}
FOLDER['table unagitated'] = 'Stability Measurements/2016-08-15 Stability Test Unagitated/Data2/'
FOLDER['table agitated'] = 'Stability Measurements/2016-08-16 Stability Test Agitated/Data/'
FOLDER['bench unagitated'] = 'Stability Measurements/2016-07-22 Stability Test/data_unagitated2/'
FOLDER['bench agitated'] = 'Stability Measurements/2016-07-22 Stability Test/data_agitated/'

def plot_stability(info_dict, title, error=0.0):
    plt.figure()
    plt.subplot(211)    
    plt.title(title)
    for test in TESTS:
        plt.errorbar(info_dict[test]['time'], info_dict[test]['r0'], label=test, yerr=error)
    plt.xlabel('Time [s]')
    plt.ylabel('Center Shift [um]')
    plt.legend(loc='best')
    plt.subplot(212)
    for test in TESTS:
        plt.errorbar(info_dict[test]['time'], info_dict[test]['diameter'], label=test, yerr=2*error)
    plt.xlabel('Time [s]')
    plt.ylabel('Diameter Shift [um]')
    plt.legend(loc='best')

def plot_diameter_stability(info_dict, title, error=0.0):
    plt.figure()
    plt.title(title)
    for test in TESTS:
        plt.errorbar(info_dict[test]['time'], info_dict[test]['diameter'], label=test, yerr=2*error)
    plt.xlabel('Time [s]')
    plt.ylabel('Diameter Shift [um]')
    plt.legend(loc='best')

def plot_center_stability(info_dict, title, error=0.0):
    plt.figure()
    plt.subplot(311)    
    plt.title(title)
    for test in TESTS:
        plt.errorbar(info_dict[test]['time'], info_dict[test]['x0'], label=test, yerr=error)
    plt.xlabel('Time [s]')
    plt.ylabel('Position Shift [um]')
    plt.legend(title='Center X', loc='best')
    plt.subplot(312)
    for test in TESTS:
        plt.errorbar(info_dict[test]['time'], info_dict[test]['y0'], label=test, yerr=error)
    plt.xlabel('Time [s]')
    plt.ylabel('Position Shift [um]')
    plt.legend(title='Center Y', loc='best')
    plt.subplot(313)
    for test in TESTS:
        plt.errorbar(info_dict[test]['time'], info_dict[test]['r0'], label=test, yerr=error)
    plt.xlabel('Time [s]')
    plt.ylabel('Position Shift [um]')
    plt.legend(title='Center R', loc='best')

if __name__ == "__main__":
    base_dict = {'x0': [], 'y0': [], 'r0': [], 'diameter': [], 'time': []}
    data_dict = {}
    for camera in CAMERAS:
        data_dict[camera] = {}


    for camera in CAMERAS:
        if camera == 'in':
            name = 'Fiber Input'
            method = 'radius'
        if camera == 'nf':
            name = 'Near Field'
            method = 'radius'
        if camera == 'ff':
            name = 'Far Field'
            method = 'gaussian'

        print name

        data_dict[camera] = {}
        for test in TESTS:
            data_dict[camera][test] = deepcopy(base_dict)

            for i in xrange(NUM_IMAGES):
                data_FOLDER = FOLDER[test] + camera + '_' + str(i).zfill(3) + '_data.p'
                obj = ImageAnalysis(image_input=None, image_data=data_FOLDER)
                y0, x0, diameter = obj.get_fiber_data(method=method, units='microns')
                if camera == 'in':
                    gauss_y, gauss_x = obj.get_fiber_center(method='gaussian', units='microns')
                    data_dict[camera][test]['x0'].append(gauss_x-x0)
                    data_dict[camera][test]['y0'].append(gauss_y-y0)
                else:
                    data_dict[camera][test]['x0'].append(x0)
                    data_dict[camera][test]['y0'].append(y0)
                data_dict[camera][test]['diameter'].append(diameter)
                data_dict[camera][test]['time'].append(obj.get_image_info('date_time'))

            for prop in ['x0', 'y0', 'diameter', 'time']:
                data_dict[camera][test][prop] = (np.array(data_dict[camera][test][prop])
                                                 - data_dict[camera][test][prop][0])
            data_dict[camera][test]['r0'] = np.sqrt(data_dict[camera][test]['x0']**2
                                                    + data_dict[camera][test]['y0']**2)
            for i, time in enumerate(data_dict[camera][test]['time']):
                data_dict[camera][test]['time'][i] = time.total_seconds()

            print test
            print 'r0 STDEV:', data_dict[camera][test]['r0'].std()
            print 'diam STDEV:', data_dict[camera][test]['diameter'].std()

        if method == 'radius':
            error = 0.1
        elif method == 'gaussian':
            error = 0.1

        error = convert_pixels_to_microns(error, pixel_size=obj.get_pixel_size(), magnification=obj.get_magnification())

        plot_stability(data_dict[camera], name + ' Stability', error)
        plt.savefig(SAVE_FOLDER + name + ' Stability.png')
        plot_center_stability(data_dict[camera], name + ' Center Stability', error)
        plt.savefig(SAVE_FOLDER + name + ' Center Stability.png')

        print

        
