from FiberProperties import ImageAnalysis, Calibration
from FiberProperties.NumpyArrayHandler import circleArray, gaussianArray, showImageArray, plotCrossSections
import numpy as np
import matplotlib.pyplot as plt

def testMethod(method):
    tol = 1
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

        mesh_grid = np.meshgrid(np.arange(image_width), np.arange(image_height))
        dark_image = [dark_level * (np.random.rand(image_height, image_width) + 0.5) for i in xrange(num_darks)]
        # noise = [0 for i in xrange(num_images)]
        noise = [np.random.rand(image_height, image_width) - 0.5 for i in xrange(num_images)]

        if method == 'edge' or method == 'radius':
            test_array = amp * circleArray(mesh_grid, x0, y0, radius, res=10) + dark_level
            test_array = [test_array + np.sqrt(test_array) * noise[i] * 30 for i in xrange(num_images)]
        elif method == 'gaussian':
            test_array = gaussianArray(mesh_grid, x0, y0, radius, amp, dark_level).reshape(image_height, image_width)
            test_array = [test_array + np.sqrt(test_array) * noise[i] * 30 for i in xrange(num_images)]

        im_obj = ImageAnalysis(test_array, Calibration(dark=dark_image), magnification=1, pixel_size=1)

        center_y, center_x, diameter = im_obj.getFiberData(method=method, units='pixels',
                                                           tol=tol, test_range=test_range)

        y_err.append(center_y - y0)
        x_err.append(center_x - x0)
        d_err.append(diameter - radius*2)

    return y_err, x_err, d_err

if __name__ == '__main__':

#=============================================================================#

    print 'Gaussian Method Testing'
    y_err, x_err, d_err = testMethod('gaussian')

    y_err = np.array(y_err)
    print 'y0 error mean:', y_err.mean(), 'y0 std:', y_err.std()
    x_err = np.array(x_err)
    print 'x0 error mean:', x_err.mean(), 'x0 std:', x_err.std()
    d_err = np.array(d_err)
    print 'diam error mean:', d_err.mean(), 'diam std:', d_err.std()

    plt.figure(1)
    plt.title('Image Analysis Accuracy')
    plt.plot(y_err, label='y0 error')
    plt.plot(x_err, label='x0 error')
    plt.plot(d_err, label='diameter error')
    plt.legend()
    plt.xlabel('Test number')
    plt.ylabel('Analysis - Actual [pixels]')
    plt.show()

#=============================================================================#

    print 'Edge Method Testing'
    y_err, x_err, d_err = testMethod('edge')

    y_err = np.array(y_err)
    print 'y0 error mean:', y_err.mean(), 'y0 std:', y_err.std()
    x_err = np.array(x_err)
    print 'x0 error mean:', x_err.mean(), 'x0 std:', x_err.std()
    d_err = np.array(d_err)
    print 'diam error mean:', d_err.mean(), 'diam std:', d_err.std()

    plt.figure(1)
    plt.title('Image Analysis Accuracy')
    plt.plot(y_err, label='y0 error')
    plt.plot(x_err, label='x0 error')
    plt.plot(d_err, label='diameter error')
    plt.legend()
    plt.xlabel('Test number')
    plt.ylabel('Analysis - Actual [pixels]')
    plt.show()

#=============================================================================#

    print 'Radius Method Testing'
    y_err, x_err, d_err = testMethod('radius')

    y_err = np.array(y_err)
    print 'y0 error mean:', y_err.mean(), 'y0 std:', y_err.std()
    x_err = np.array(x_err)
    print 'x0 error mean:', x_err.mean(), 'x0 std:', x_err.std()
    d_err = np.array(d_err)
    print 'diam error mean:', d_err.mean(), 'diam std:', d_err.std()

    plt.figure(1)
    plt.title('Image Analysis Accuracy')
    plt.plot(y_err, label='y0 error')
    plt.plot(x_err, label='x0 error')
    plt.plot(d_err, label='diameter error')
    plt.legend()
    plt.xlabel('Test number')
    plt.ylabel('Analysis - Actual [pixels]')
    plt.show()