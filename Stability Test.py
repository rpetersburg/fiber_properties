import matplotlib.pyplot as plt
import numpy as np
from ImageAnalysis import ImageAnalysis
from Calibration import Calibration
import re

import os
from multiprocessing import Pool
from functools import partial

def saveInputData(image, in_calibration):
    obj = ImageAnalysis(image, in_calibration)
    obj.setFiberData(method='edge')
    obj.setFiberData(method='radius', tol=0.1, test_range=5)
    obj.setFiberData(method='gaussian')
    obj.setFiberCentroid(method='full')
    file_name = re.split('/|\.', image)[-2]
    obj.saveData(file_name=file_name)
    print file_name + ' saved'

def saveNearFieldData(image, nf_calibration):
    obj = ImageAnalysis(image, nf_calibration)
    obj.setFiberData(method='edge')
    obj.setFiberData(method='radius', tol=0.1, test_range=5)
    obj.setFiberCentroid(method='full')
    file_name = re.split('/|\.', image)[-2]
    obj.saveData(file_name=file_name)
    print file_name + ' saved'

def saveFarFieldData(image, ff_calibration):
    obj = ImageAnalysis(image, ff_calibration, magnification=1)
    obj.setFiberData(method='gaussian')
    obj.setFiberCentroid(method='full')
    file_name = re.split('/|\.', image)[-2]
    obj.saveData(file_name=file_name)
    print file_name + ' saved'

if __name__ == '__main__':
    base_folder = 'Stability Measurements/2016-08-15 Stability Test Unagitated/'
    ambient_folder = base_folder + 'Ambient/'
    dark_folder = base_folder + 'Dark/'
    flat_folder = base_folder + 'Flat/'
    folder = base_folder + 'Images/'
    ext = '.fit'

    in_calibration = Calibration(dark=[dark_folder + 'in_' + str(i).zfill(3) + ext for i in xrange(10)],
                                 ambient=[ambient_folder + 'in_' + str(i).zfill(3) + ext for i in xrange(10)],
                                 flat=[flat_folder + 'in_' + str(i).zfill(3) + ext for i in xrange(8)])
    print 'IN calibration initialized'
    nf_calibration = Calibration(dark=[dark_folder + 'nf_' + str(i).zfill(3) + ext for i in xrange(10)],
                                 ambient=[ambient_folder + 'nf_' + str(i).zfill(3) + ext for i in xrange(10)],
                                 flat=[flat_folder + 'nf_' + str(i).zfill(3) + ext for i in xrange(8)])
    print 'NF calibration initialized'
    ff_calibration = Calibration(dark=[dark_folder + 'ff_' + str(i).zfill(3) + ext for i in xrange(10)],
                                 ambient=[ambient_folder + 'ff_' + str(i).zfill(3) + ext for i in xrange(10)])
    print 'FF calibration initialized'

    num_images = 100
    in_images = [folder + 'in_' + str(i).zfill(3) + ext for i in xrange(num_images)]
    nf_images = [folder + 'nf_' + str(i).zfill(3) + ext for i in xrange(num_images)]
    ff_images = [folder + 'ff_' + str(i).zfill(3) + ext for i in xrange(num_images)]

    pool = Pool(processes=2)
    func = partial(saveInputData, in_calibration=in_calibration)
    pool.map(func, in_images)
    print 'IN data initialized'
    func = partial(saveNearFieldData, nf_calibration=nf_calibration)
    pool.map(func, nf_images)
    print 'NF data initialized'
    func = partial(saveFarFieldData, ff_calibration=ff_calibration)
    pool.map(func, ff_images)
    print 'FF data initialized'

    # for image in in_images:
    #     saveInputData(image)
    # print 'IN data initialized'

    # for image in nf_images:
    #     saveNearFieldData(image)
    # print 'NF data initialized'

    # for image in ff_images:
    #     saveFarFieldData(image)
    # print 'FF data initialized'

    # plt.figure(1)
    # plt.title('Near Field Diameter Stability')
    # plt.plot(nf['agitated']['diameter'], label='agitated')
    # plt.plot(nf['unagitated']['diameter'], label='unagitated')
    # plt.xlabel('Image number')
    # plt.ylabel('Diameter [um]')
    # plt.legend()

    # plt.figure(2)
    # plt.subplot(211)
    # plt.plot(nf['agitated']['center_x'], label='agitated')
    # plt.plot(nf['unagitated']['center_x'], label='unagitated')
    # plt.xlabel('Image number')
    # plt.ylabel('Position [um]')
    # plt.legend(title='Center X')
    # plt.title('Near Field Center Stability')
    # plt.subplot(212)
    # plt.plot(nf['agitated']['center_y'], label='agitated')
    # plt.plot(nf['unagitated']['center_y'], label='unagitated')
    # plt.xlabel('Image number')
    # plt.ylabel('Position [um]')
    # plt.legend(title='Center Y')

    # plt.figure(3)
    # plt.subplot(211)
    # plt.plot(nf['agitated']['centroid_x'], label='agitated')
    # plt.plot(nf['unagitated']['centroid_x'], label='unagitated')
    # plt.xlabel('Image number')
    # plt.ylabel('Position [um]')
    # plt.legend(title='Centroid X')
    # plt.title('Near Field Centroid Stability')
    # plt.subplot(212)
    # plt.plot(nf['agitated']['centroid_y'], label='agitated')
    # plt.plot(nf['unagitated']['centroid_y'], label='unagitated')
    # plt.xlabel('Image number')
    # plt.ylabel('Position [um]')
    # plt.legend(title='Centroid Y')

    # plt.figure(4)
    # plt.title('Far Field Diameter Stability')
    # plt.plot(ff['agitated']['diameter'], label='agitated')
    # plt.plot(ff['unagitated']['diameter'], label='unagitated')
    # plt.xlabel('Image number')
    # plt.ylabel('Diameter [um]')
    # plt.legend()

    # plt.figure(5)
    # plt.subplot(211)
    # plt.plot(ff['agitated']['center_x'], label='agitated')
    # plt.plot(ff['unagitated']['center_x'], label='unagitated')
    # plt.xlabel('Image number')
    # plt.ylabel('Position [um]')
    # plt.legend(title='Center X [um]')
    # plt.title('Far Field Center Stability')
    # plt.subplot(212)
    # plt.plot(ff['agitated']['center_y'], label='agitated')
    # plt.plot(ff['unagitated']['center_y'], label='unagitated')
    # plt.xlabel('Image number')
    # plt.ylabel('Position [um]')
    # plt.legend(title='Center Y [um]')

    # plt.figure(6)
    # plt.subplot(211)
    # plt.plot(ff['agitated']['centroid_x'], label='agitated')
    # plt.plot(ff['unagitated']['centroid_x'], label='unagitated')
    # plt.xlabel('Image number')
    # plt.ylabel('Position [um]')
    # plt.legend(title='Centroid X [um]')
    # plt.title('Far Field Centroid Stability')
    # plt.subplot(212)
    # plt.plot(ff['agitated']['centroid_y'], label='agitated')
    # plt.plot(ff['unagitated']['centroid_y'], label='unagitated')
    # plt.xlabel('Image number')
    # plt.ylabel('Position [um]')
    # plt.legend(title='Centroid Y [um]')

    # plt.show()
