import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy
from ImageAnalysis import ImageAnalysis
from Calibration import Calibration

base_folder = '2016-07-22/'
ambient_folder = base_folder + 'ambient/'
dark_folder = base_folder + 'dark/'
flat_folder = base_folder + 'flat/'
agitated_folder = base_folder + 'stability_agitated/'
unagitated_folder = base_folder + 'stability_unagitated/'
file_extension = '.fit'

nf = {}
ff = {}

nf['calibration'] = Calibration([dark_folder + 'nf_dark_' + str(i).zfill(3) + file_extension for i in xrange(10)],
                                [flat_folder + 'nf_flat_' + str(i) + '_1ms' + file_extension for i in xrange(8)],
                                [ambient_folder + 'nf_ambient_' + str(i).zfill(3) + '_0.001' + file_extension for i in xrange(10)])
print 'NF calibration initialized'
ff['calibration'] = Calibration([dark_folder + 'ff_dark_' + str(i).zfill(3) + file_extension for i in xrange(10)],
                                None,
                                [ambient_folder + 'ff_ambient_' + str(i).zfill(3) + '_0.001' + file_extension for i in xrange(10)])
print 'FF calibration initialized'

empty_data = {'images': [], 'diameter': [], 'center_x': [], 'center_y': [], 'centroid_x': [], 'centroid_y': []}
nf['agitated'] = deepcopy(empty_data)
nf['unagitated'] = deepcopy(empty_data)
ff['agitated'] = deepcopy(empty_data)
ff['unagitated'] = deepcopy(empty_data)

image_range = xrange(100)
nf['agitated']['images'] = [agitated_folder + 'nf_stability_' + str(i).zfill(3) + '_0.002' + file_extension for i in image_range]
nf['unagitated']['images'] = [unagitated_folder + 'nf_stability_' + str(i).zfill(3) + '_0.002' + file_extension for i in image_range]

ff['agitated']['images'] = [agitated_folder + 'ff_stability,fit_' + str(i).zfill(3) + '_0.001' + file_extension for i in image_range]
ff['unagitated']['images'] = [unagitated_folder + 'ff_stability,fit_' + str(i).zfill(3) + '_0.001' + file_extension for i in image_range]

for test in ['agitated', 'unagitated']:
    for image in nf[test]['images']:
        obj = ImageAnalysis(image, nf['calibration'])
        nf[test]['diameter'].append(obj.getFiberDiameter(method='edge', units='microns'))
        y, x = obj.getFiberCenter(method='edge', show_image=False)
        nf[test]['center_x'].append(x)
        nf[test]['center_y'].append(y)
        y, x = obj.getFiberCentroid()
        nf[test]['centroid_x'].append(x)
        nf[test]['centroid_y'].append(y)
        print image

write_file = open('2016-07-25/stability_data.txt', 'w')
write_file.write(str(nf))
write_file.close()

# for test in ['agitated', 'unagitated']:
#     for image in ff[test]['images']:
#         obj = ImageAnalysis(image, f['calibration'])
#         ff[test]['diameter'].append(obj.getFiberDiameter(method='edge', units='microns'))
#         y, x = obj.getFiberCenter(method='edge', show_image=False)
#         ff[test]['center_x'].append(x)
#         ff[test]['center_y'].append(y)
#         y, x = obj.getFiberCentroid()
#         ff[test]['centroid_x'].append(x)
#         ff[test]['centroid_y'].append(y)
#         print image
#     print 'NF', test, 'objects initialized'

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

# plt.figure(4)
# plt.plot(ff_agitated_diameters, label='agitated')
# plt.plot(ff_unagitated_diameters, label='unagitated')
# plt.title('Far Field Diameter Stability')
# plt.xlabel('Image number')
# plt.ylabel('Diameter [um]')
# plt.legend()

# plt.figure(5)
# plt.plot(ff_agitated_centers, label='agitated')
# plt.plot(ff_unagitated_centers, label='unagitated')
# plt.title('Far Field Center Stability')
# plt.xlabel('Image number')
# plt.ylabel('Center radial distance [um]')
# plt.legend()

# plt.figure(6)
# plt.plot(ff_agitated_centroids, label='agitated')
# plt.plot(ff_unagitated_centroids, label='unagitated')
# plt.title('Far Field Centroid Stability')
# plt.xlabel('Image number')
# plt.ylabel('Centroid radial distance [um]')
# plt.legend()

plt.show()
