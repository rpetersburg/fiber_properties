import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy
from ImageAnalysis import ImageAnalysis
from Calibration import Calibration

base_folder = '2016-07-26/single/100um/'
ambient_folder = '2016-07-26/ambient/100um/'
dark_folder = '2016-07-26/dark/'
agitated_folder = base_folder + 'agitated/'
unagitated_folder = base_folder + 'unagitated/'
ext = '.fit'

nf = {}
ff = {}

nf['calibration'] = Calibration([dark_folder + 'nf_dark_' + str(i).zfill(3) + ext for i in xrange(10)],
                                None,
                                [ambient_folder + 'nf_ambient_' + str(i).zfill(3) + '_0.1' + ext for i in xrange(10)])
print 'NF calibration initialized'
ff['calibration'] = Calibration([dark_folder + 'ff_dark_' + str(i).zfill(3) + ext for i in xrange(10)],
                                None,
                                [ambient_folder + 'ff_ambient_' + str(i).zfill(3) + '_0.1' + ext for i in xrange(10)])
print 'FF calibration initialized'

empty_data = {'images': [], 'diameter': [], 'center_x': [], 'center_y': [], 'centroid_x': [], 'centroid_y': []}
nf['agitated'] = deepcopy(empty_data)
nf['unagitated'] = deepcopy(empty_data)
ff['agitated'] = deepcopy(empty_data)
ff['unagitated'] = deepcopy(empty_data)

image_range = xrange(20)
nf['agitated']['images'] = [agitated_folder + 'nf_agitated_' + str(i).zfill(3) + ext for i in image_range]
nf['unagitated']['images'] = [unagitated_folder + 'nf_unagitated_' + str(i).zfill(3) + ext for i in image_range]

ff['agitated']['images'] = [agitated_folder + 'ff_agitated_' + str(i).zfill(3) + ext for i in image_range]
ff['unagitated']['images'] = [unagitated_folder + 'ff_unagitated_' + str(i).zfill(3) + ext for i in image_range]

for test in ['agitated', 'unagitated']:
    for image in nf[test]['images']:
        obj = ImageAnalysis(image, nf['calibration'])
        y, x, diameter = obj.getFiberData(method='radius', tol=0.1, test_range=5, units='microns', show_image=True)
        nf[test]['diameter'].append(diameter)
        nf[test]['center_x'].append(x)
        nf[test]['center_y'].append(y)
        y, x = obj.getFiberCentroid()
        nf[test]['centroid_x'].append(x)
        nf[test]['centroid_y'].append(y)
        print image
    print 'NF', test, 'objects initialized'

for test in ['agitated', 'unagitated']:
    for image in ff[test]['images']:
        obj = ImageAnalysis(image, ff['calibration'])
        ff[test]['diameter'].append(obj.getFiberDiameter(method='gaussian', units='microns'))
        y, x = obj.getFiberCenter(method='gaussian', show_image=False, units='microns')
        ff[test]['center_x'].append(x)
        ff[test]['center_y'].append(y)
        y, x = obj.getFiberCentroid()
        ff[test]['centroid_x'].append(x)
        ff[test]['centroid_y'].append(y)
        print image
    print 'FF', test, 'objects initialized'

write_file = open('2016-07-26/stability_data.txt', 'w')
write_file.write(str(ff))
write_file.close()

plt.figure(1)
plt.title('Near Field Diameter Stability')
plt.plot(nf['agitated']['diameter'], label='agitated')
plt.plot(nf['unagitated']['diameter'], label='unagitated')
plt.xlabel('Image number')
plt.ylabel('Diameter [um]')
plt.legend()

plt.figure(2)
plt.subplot(211)
plt.plot(nf['agitated']['center_x'], label='agitated')
plt.plot(nf['unagitated']['center_x'], label='unagitated')
plt.xlabel('Image number')
plt.ylabel('Position [um]')
plt.legend(title='Center X')
plt.title('Near Field Center Stability')
plt.subplot(212)
plt.plot(nf['agitated']['center_y'], label='agitated')
plt.plot(nf['unagitated']['center_y'], label='unagitated')
plt.xlabel('Image number')
plt.ylabel('Position [um]')
plt.legend(title='Center Y')

plt.figure(3)
plt.subplot(211)
plt.plot(nf['agitated']['centroid_x'], label='agitated')
plt.plot(nf['unagitated']['centroid_x'], label='unagitated')
plt.xlabel('Image number')
plt.ylabel('Position [um]')
plt.legend(title='Centroid X')
plt.title('Near Field Centroid Stability')
plt.subplot(212)
plt.plot(nf['agitated']['centroid_y'], label='agitated')
plt.plot(nf['unagitated']['centroid_y'], label='unagitated')
plt.xlabel('Image number')
plt.ylabel('Position [um]')
plt.legend(title='Centroid Y')

plt.figure(4)
plt.title('Far Field Diameter Stability')
plt.plot(ff['agitated']['diameter'], label='agitated')
plt.plot(ff['unagitated']['diameter'], label='unagitated')
plt.xlabel('Image number')
plt.ylabel('Diameter [um]')
plt.legend()

plt.figure(5)
plt.subplot(211)
plt.plot(ff['agitated']['center_x'], label='agitated')
plt.plot(ff['unagitated']['center_x'], label='unagitated')
plt.xlabel('Image number')
plt.ylabel('Position [um]')
plt.legend(title='Center X [um]')
plt.title('Far Field Center Stability')
plt.subplot(212)
plt.plot(ff['agitated']['center_y'], label='agitated')
plt.plot(ff['unagitated']['center_y'], label='unagitated')
plt.xlabel('Image number')
plt.ylabel('Position [um]')
plt.legend(title='Center Y [um]')

plt.figure(6)
plt.subplot(211)
plt.plot(ff['agitated']['centroid_x'], label='agitated')
plt.plot(ff['unagitated']['centroid_x'], label='unagitated')
plt.xlabel('Image number')
plt.ylabel('Position [um]')
plt.legend(title='Centroid X [um]')
plt.title('Far Field Centroid Stability')
plt.subplot(212)
plt.plot(ff['agitated']['centroid_y'], label='agitated')
plt.plot(ff['unagitated']['centroid_y'], label='unagitated')
plt.xlabel('Image number')
plt.ylabel('Position [um]')
plt.legend(title='Centroid Y [um]')

plt.show()
