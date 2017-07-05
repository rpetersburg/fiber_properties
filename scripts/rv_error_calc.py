from fiber_properties import FiberImage
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

NUM_IMAGES = 1
CASE = 3
FOLDER = '/Users/Dominic/Box Sync/Fiber_Characterization/Image Analysis/data/modal_noise/rv_error/'
METHOD = 'radius'
CAMERA = ['nf', 'ff']

if CASE == 1:
    FOLDER += 'coupled_agitation/radius_pkls/'
if CASE == 2:
    FOLDER += 'LED/radius_pkls/'
if CASE == 3:
    FOLDER += 'slow_agitation/radius_pkls/'


def find_rv_error(folder=FOLDER, num_images=NUM_IMAGES, camera=CAMERA, meth=METHOD):

    for cam in camera:
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
        for i in enumerate(center_x):
            i = i[0]  # I is a tuple, [0] is the index
            calc = np.sqrt(center_x[i]**2 + center_y[i]**2) * (np.cos(np.arctan(center_y[i]/center_x[i]) - np.pi/3))
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

        # Plot #
        center_line = plt.plot(center, color='g', label='$\sigma_{cen}=%.2f$\n$\sigma_{rv}=%.2f$' % (center_std, rv_std_all))
        avg_line = plt.plot(xrange(5, 300, 10), center_avg, color='r', label='$\sigma_{avg}=%.2f$\n$\sigma_{rv}=%.2f$' % (center_avg_std, rv_std_avg))

        plt.ylabel('Center drift ($\mu m$)')
        plt.xlabel('Frame number')

        plt.legend(loc='lower right')

        # Save #
        plt.savefig(folder + 'plots/rv_error_plots/%s_%s_rv_error.png' % (cam, meth), bbox_inches='tight')
        print('Saved figure to %splots/rv_error_plots/' % str(folder))
        plt.close()

if __name__ == '__main__':
    find_rv_error()
