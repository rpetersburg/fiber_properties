from fiber_properties import FiberImage
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

NUM_IMAGES = 1
CASE = 3
FOLDER = '/Users/Dominic/Box Sync/Fiber_Characterization/Image Analysis/data/modal_noise/rv_error/'
METHOD = 'full'
CAMERAS = ['ff']

if CASE == 1:
    FOLDER += 'coupled_agitation/'
if CASE == 2:
    FOLDER += 'LED/'
if CASE == 3:
    FOLDER += 'slow_agitation/'


def find_rv_error(folder=FOLDER, num_images=NUM_IMAGES, camera=CAMERAS, meth=METHOD):

    for cam in camera:
        print('Saving to Folder: %s' % folder)
        print('Camera: %s' % str(cam))
        print('Method: %s' % str(meth))
        center_x = []
        center_y = []
        for i in xrange(0, 300, num_images):
            im_obj = FiberImage(folder + cam + '_' + str(i).zfill(3) + '_obj.pkl')
            c = im_obj.get_fiber_center(method=meth, units='microns') - im_obj.get_fiber_centroid(method=meth, units='microns')
            center_x.append(c.x)
            center_y.append(c.y)
            print('Getting center for image set %s...' % i)
            print(c)

        cen_x = np.array(center_x)
        cen_y = np.array(center_y)
        # Compute drift #
        x_median = np.median(center_x)
        y_median = np.median(center_y)
        center_x = [x - x_median for x in cen_x]
        center_y = [y - y_median for y in cen_y]

        center = []
        for c_x, c_y in zip(center_x, center_y):
            calc = np.sqrt(c_x**2 + c_y**2) * np.cos(np.arctan(c_y/c_x) + (np.pi/6) + (np.pi/2)*(1 - np.sign(c_x)))
            center.append(calc)

        center_std = np.std(center)

        if cam is 'nf':
            diameter = 100
        if cam is 'ff':
            print('Getting ff diameter')
            diameter = im_obj.get_fiber_diameter(method=meth, units='microns')

        # Make average line #
        center_avg = []
        itr = range(10)
        for i in xrange(0, 300, 10):
            center_to_avg = []
            for num in itr:
                center_to_avg.append(center[num])
            center_avg.append(np.average(center_to_avg))
            itr = [x + 10 for x in itr]

        center_avg_std = np.std(center_avg)

        # Compute RV error #
        rv_std_all = (3 * 10**8) * (center_std) / ((150000)*(diameter))
        rv_std_avg = (3 * 10**8) * (center_avg_std) / ((150000)*(diameter))

        # Convert to m/s #
        center_ms = [(3 * 10**8) * (x) / ((150000)*(diameter)) for x in center]
        center_avg_ms = [(3 * 10**8) * (x) / ((150000)*(diameter)) for x in center_avg]

        # Plot #
        center_line = plt.plot(center_ms, color='g', label='$\sigma_{rv}=%.2f$' % (rv_std_all))
        avg_line = plt.plot(xrange(5, 300, 10), center_avg_ms, color='r', label='$\sigma_{rv}=%.2f$' % (rv_std_avg))

        plt.ylabel('Center drift (m/s)')
        plt.xlabel('Frame number')

        plt.legend(loc='lower right')

        # Save #
        plt.savefig(folder + 'plots/rv_error_plots/%s_%s_rv_error.png' % (cam, meth), bbox_inches='tight')
        print('Saved figure to %splots/rv_error_plots/' % str(folder))
        plt.close()

if __name__ == '__main__':
    find_rv_error()
