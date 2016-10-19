import re
import os
from ImageAnalysis import ImageAnalysis

def saveInputData(image, in_calibration): 
    file_name = 'Data/' + re.split('/|\.', image)[-2]

    if file_name + '_data.txt' in os.listdir('/'.join(image.split('/')[:-1])):
        print file_name + ' already saved'
        return

    obj = ImageAnalysis(image, in_calibration, camera='in')
    print file_name + ' initialized'

    obj.setFiberData(method='edge')
    obj.setFiberData(method='radius', tol=0.25, test_range=5)
    obj.setFiberData(method='gaussian')
    obj.setFiberCentroid(method='full')

    obj.saveData(file_name=file_name)
    print file_name + ' saved, diameter:', obj.getFiberDiameter(method='radius', units='microns')

def saveNearFieldData(image, nf_calibration):
    file_name = 'Data/' + re.split('/|\.', image)[-2]

    if file_name + '_data.txt' in os.listdir('/'.join(image.split('/')[:-1])):
        print file_name + ' already saved'
        return

    obj = ImageAnalysis(image, nf_calibration, camera='nf', threshold=500)
    print file_name + ' initialized'

    obj.setFiberData(method='edge')
    obj.setFiberData(method='radius', tol=0.25, test_range=5)
    obj.setFiberCentroid(method='full')

    obj.saveData(file_name=file_name)
    print file_name + ' saved, diameter:', obj.getFiberDiameter(method='radius', units='microns')

def saveFarFieldData(image, ff_calibration):
    file_name = 'Data/' + re.split('/|\.', image)[-2]

    if file_name + '_data.txt' in os.listdir('/'.join(image.split('/')[:-1])):
        print file_name + ' already saved'
        return

    obj = ImageAnalysis(image, ff_calibration, camera='ff', magnification=1)
    print file_name + ' initialized'

    obj.setFiberData(method='gaussian')
    obj.setFiberCentroid(method='full')

    obj.saveData(file_name=file_name)
    print file_name + ' saved, diameter:', obj.getFiberDiameter(method='gaussian', units='microns')

if __name__ == '__main__':
    from Calibration import Calibration

    PARALLELIZE = True
    
    base_folder = 'Stability Measurements/2016-07-22 Stability Test/'
    ambient_folder = base_folder + 'Ambient/'
    dark_folder = base_folder + 'Dark/'
    flat_folder = base_folder + 'Flat/'
    folder = base_folder + 'stability_unagitated/'
    ext = '.fit'

    in_calibration = Calibration(dark=[dark_folder + 'in_' + str(i).zfill(3) + ext for i in xrange(10)],
                                 ambient=[ambient_folder + 'in_' + str(i).zfill(3) + ext for i in xrange(10)],
                                 flat=[flat_folder + 'in_' + str(i).zfill(3) + ext for i in xrange(8)])
    print 'IN calibration initialized'
    nf_calibration = Calibration(dark=[dark_folder + 'nf_' + str(i).zfill(3) + ext for i in xrange(10)],
                                 ambient=[ambient_folder + 'nf_' + str(i).zfill(3) + '_0.001' + ext for i in xrange(10)],
                                 flat=[flat_folder + 'nf_' + str(i).zfill(3) + ext for i in xrange(8)])
    print 'NF calibration initialized'
    ff_calibration = Calibration(dark=[dark_folder + 'ff_' + str(i).zfill(3) + ext for i in xrange(10)],
                                 ambient=[ambient_folder + 'ff_' + str(i).zfill(3) + '_0.001' + ext for i in xrange(10)])
    print 'FF calibration initialized'

    num_images = 100
    in_images = [folder + 'in_' + str(i).zfill(3) + ext for i in xrange(num_images)]
    nf_images = [folder + 'nf_' + str(i).zfill(3) + ext for i in xrange(num_images)]
    ff_images = [folder + 'ff_' + str(i).zfill(3) + ext for i in xrange(num_images)]

    if PARALLELIZE:
        from multiprocessing import Pool
        from functools import partial
        pool = Pool(processes=6)
        pool.map(partial(saveInputData, in_calibration=in_calibration), in_images)
        print 'IN data initialized'
        pool.map(partial(saveNearFieldData, nf_calibration=nf_calibration), nf_images)
        print 'NF data initialized'
        pool.map(partial(saveFarFieldData, ff_calibration=ff_calibration), ff_images)
        print 'FF data initialized'

    else:    
        for image in in_images:
            saveInputData(image, in_calibration)
        print 'IN data initialized'

        for image in nf_images:
            saveNearFieldData(image, nf_calibration)
        print 'NF data initialized'

        for image in ff_images:
            saveFarFieldData(image, ff_calibration)
        print 'FF data initialized'