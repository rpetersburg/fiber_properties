from ImageAnalysis import ImageAnalysis
from Calibration import Calibration
from NumpyArrayHandler import NumpyArrayHandler as NAH
import numpy as np
import matplotlib.pyplot as plt

tol = 0.1
test_range = 1
factor = 1.0

y_err = []
x_err = []
d_err = []
    
for i in xrange(10):
    image_height = 2000
    image_width = 2000
    x0 = image_width * 0.5 + np.random.rand()
    y0 = image_height * 0.5 + np.random.rand()
    radius = min(image_height, image_width) * 0.25 + np.random.rand()
    amp = 10000
    dark_level = 500
    num_images = 1
    num_darks = 10

    dark_image = [dark_level * (np.random.rand(image_height, image_width) + 0.5) for i in xrange(num_darks)]

    mesh_grid = np.meshgrid(np.arange(image_width), np.arange(image_height))
    noise = [np.random.rand(image_height, image_width) + 0.5 for i in xrange(num_images)]

    circle_array = [(amp * NAH.circleArray(mesh_grid, x0, y0, radius, res=10) + dark_level) * noise[i] for i in xrange(num_images)]

    im_obj = ImageAnalysis(circle_array, Calibration(dark=dark_image), pixel_size=1,
                           camera='nf', threshold=None, kernel_size=1)

    center_y, center_x, diameter = im_obj.getFiberData(method='radius', tol=tol,
                                                       test_range=test_range,
                                                       show_image=False, units='pixels')
    y_err.append(center_y - y0)
    x_err.append(center_x - x0)
    d_err.append(diameter - radius*2)

y_err = np.array(y_err)
print 'y0 error mean:', np.abs(y_err).mean(), 'y0 std:', y_err.std()
x_err = np.array(x_err)
print 'x0 error mean:', np.abs(x_err).mean(), 'x0 std:', x_err.std()
d_err = np.array(d_err)
print 'diam error mean:', np.abs(d_err).mean(), 'diam std:', d_err.std()

plt.figure(1)
plt.title('Image Analysis Accuracy')
plt.plot(y_err, label='y0 error')
plt.plot(x_err, label='x0 error')
plt.plot(d_err, label='diameter error')
plt.legend()
plt.xlabel('Test number')
plt.ylabel('Analysis - Actual [microns]')
plt.show()



# print 'Actual Data:'
# print 'Diameter:', radius * 2.0, 'pixels'
# print 'Center Row:', y0, 'Center Column:', x0
# print
# print 'Analyzed Data:'
# print 'Centroid:'
# centroid_row, centroid_column = im_obj.getFiberCentroid(factor)
# print 'Centroid Row:', centroid_row, 'Centroid Column:', centroid_column
# print
# print 'Edge:'
# center_y, center_x = im_obj.getFiberCenter(method='edge', show_image=True)
# print 'Diameter:', im_obj.getFiberDiameter(method='edge', show_image=False, units='pixels'), 'pixels'
# print 'Center Row:', center_y, 'Center Column:', center_x
# print
# print 'Radius:'
# center_y, center_x = im_obj.getFiberCenter(method='radius', tol=tol, show_image=True, test_range=test_range)
# print 'Diameter:', im_obj.getFiberDiameter(method='radius',show_image=False, units='pixels'), 'pixels'
# print 'Center Row:', center_y, 'Center Column:', center_x
# print

# gauss_array = [NAH.gaussianArray(mesh_grid, x0, y0, radius, amp, dark_level).reshape(image_height, image_width) * noise[i] for i in xrange(num_images)]
# im_obj = ImageAnalysis(gauss_array, Calibration(dark=dark_image), pixel_size=1,
#                        camera='ff', threshold=None, kernel_size=1)
# im_obj.showImageArray()

# print 'Actual Data:'
# print 'Diameter:', radius * 2.0, 'pixels'
# print 'Center Row:', y0, 'Center Column:', x0
# print
# print 'Analyzed Data:'
# print 'Centroid:'
# centroid_row, centroid_column = im_obj.getFiberCentroid(factor)
# print 'Centroid Row:', centroid_row, 'Centroid Column:', centroid_column
# print
# print 'Gaussian:'
# center_y, center_x = im_obj.getFiberCenter(method='gaussian', show_image=True)
# print 'Dimaeter:', im_obj.getFiberDiameter(method='gaussian', show_image=False, units='pixels'), 'pixels'
# print 'Center Row:', center_y, 'Center Column:', center_x
# print 