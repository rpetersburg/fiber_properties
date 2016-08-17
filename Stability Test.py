import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy
from ImageAnalysis import ImageAnalysis
from Calibration import Calibration

base_folder = 'Stability Measurements/2016-08-15 Stability Test Unagitated/'
ambient_folder = base_folder + 'Ambient/'
dark_folder = base_folder + 'Dark/'
folder = base_folder + 'Images/'
ext = '.fit'

in_calibration = Calibration(dark=[dark_folder + 'in_' + str(i).zfill(3) + ext for i in xrange(10)],
                             ambient=[ambient_folder + 'in_' + str(i).zfill(3) + ext for i in xrange(10)])
print 'IN calibration initialized'
nf_calibration = Calibration(dark=[dark_folder + 'nf_' + str(i).zfill(3) + ext for i in xrange(10)],
                             ambient=[ambient_folder + 'nf_' + str(i).zfill(3) + ext for i in xrange(10)])
print 'NF calibration initialized'
ff_calibration = Calibration(dark=[dark_folder + 'ff_' + str(i).zfill(3) + ext for i in xrange(10)],
                             ambient=[ambient_folder + 'ff_' + str(i).zfill(3) + ext for i in xrange(10)])
print 'FF calibration initialized'

num_images = 10
in_images = [folder + 'in_' + str(i).zfill(3) + ext for i in xrange(num_images)]
nf_images = [folder + 'nf_' + str(i).zfill(3) + ext for i in xrange(num_images)]
ff_images = [folder + 'ff_' + str(i).zfill(3) + ext for i in xrange(num_images)]

for num, image in enumerate(in_images):
    obj = ImageAnalysis(image, in_calibration)
    obj.setFiberCenter(method='radius', tol=0.1, test_range=10, show_image=True)
    obj.setFiberCenter(method='gaussian', tol=0.1, test_range=10, show_image=True)
    obj.setFiberCentroid(method=)
    obj.setFiberData(method='gaussian', show_image=True)
    obj.saveData(file_name='in_' + str(num).zfill(3) + '_obj_data')
print 'IN data initialized'

for num, image in enumerate(nf_images):
    obj = ImageAnalysis(image, nf_calibration)
    obj.setFiberData(method='radius', tol=0.1, test_range=10, show_image=True)
    obj.saveData(file_name='nf_' + str(num).zfill(3) + '_obj_data')
print 'NF data initialized'

for num, image in enumerate(ff_images):
    obj = ImageAnalysis(image, ff_calibration, magnification=1)
    obj.setFiberData(method='gaussian', show_image=True)
    obj.saveData(file_name='ff_' + str(num).zfill(3) + 'obj_data')
print 'FF data initialized'

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
