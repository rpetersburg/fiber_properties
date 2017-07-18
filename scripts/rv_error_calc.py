from fiber_properties import FiberImage
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

PLOT_PER_10 = True
NUM_IMAGES = 1
CASE = 2
FOLDER = '/Users/Dominic/Box Sync/Fiber_Characterization/Image Analysis/data/modal_noise/rv_error/'
METHOD = 'full'
CAMERAS = ['nf', 'ff']

if CASE == 1:
    FOLDER += 'coupled_agitation/'
    ANGLE_NF = np.pi/6
    ANGLE_FF = 2*(np.pi)/3
if CASE == 2:
    FOLDER += 'LED/'
    ANGLE_NF = np.pi/6
    ANGLE_FF = 2*(np.pi)/3
if CASE == 3:
    FOLDER += 'slow_agitation/'
    ANGLE_NF = np.pi/6
    ANGLE_FF = 2*(np.pi)/3
if CASE == 4:
    FOLDER += 'coupled_ag_new/'
    ANGLE_NF = 0.611
    ANGLE_FF = 2.53


def find_rv_error(folder=FOLDER, num_images=NUM_IMAGES, camera=CAMERAS, meth=METHOD, per_10=PLOT_PER_10, angle_nf=ANGLE_NF, angle_ff=ANGLE_FF):

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
            if cam is 'nf':
                calc = np.sqrt(c_x**2 + c_y**2) * np.cos(np.arctan(c_y/c_x) + (angle_nf) + (np.pi/2)*(1 - np.sign(c_x)))
                center.append(calc)
            if cam is 'ff':
                calc = np.sqrt(c_x**2 + c_y**2) * np.cos(np.arctan(c_y/c_x) + (angle_ff) + (np.pi/2)*(1 - np.sign(c_x)))
                center.append(calc)

        # Make average line #
        if per_10:
            num = 30
            center_xavg = []
            center_yavg = []
            for i in xrange(0, 300, num):
                avg_file = FiberImage(folder + cam + '_' + str(i).zfill(3) + '_' + str(i+num-1).zfill(3) + '_obj.pkl')
                a = avg_file.get_fiber_center(method=meth, units='microns') - avg_file.get_fiber_centroid(method=meth, units='microns')
                center_xavg.append(a.x)
                center_yavg.append(a.y)
                print('Getting center for image set %s...' % i)
                print(a)

            cen_xavg = np.array(center_xavg)
            cen_yavg = np.array(center_yavg)
            # Compute drift #
            xavg_median = np.median(center_xavg)
            yavg_median = np.median(center_yavg)
            center_xavg = [x - xavg_median for x in cen_xavg]
            center_yavg = [y - yavg_median for y in cen_yavg]

            center_avg = []
            for c_x, c_y in zip(center_xavg, center_yavg):
                if cam is 'nf':
                    calc = np.sqrt(c_x**2 + c_y**2) * np.cos(np.arctan(c_y/c_x) + (angle_nf) + (np.pi/2)*(1 - np.sign(c_x)))
                    center_avg.append(calc)
                if cam is 'ff':
                    calc = np.sqrt(c_x**2 + c_y**2) * np.cos(np.arctan(c_y/c_x) + (angle_ff) + (np.pi/2)*(1 - np.sign(c_x)))
                    center_avg.append(calc)

        else:
            center_avg = []
            itr = range(10)
            for i in xrange(0, 300, 10):
                center_to_avg = []
                for num in itr:
                    center_to_avg.append(center[num])
                center_avg.append(np.average(center_to_avg))
                itr = [x + 10 for x in itr]

        # Compute std #
        center_std = np.std(center)
        center_avg_std = np.std(center_avg)

        # Get diameter #
        if cam is 'nf':
            diameter = 100
            diameter_avg = 100
        if cam is 'ff':
            print('Getting ff diameter')
            if per_10:
                diameter = im_obj.get_fiber_diameter(method=meth, units='microns')
                diameter_avg = avg_file.get_fiber_diameter(method=meth, units='microns')
            else:
                diameter = im_obj.get_fiber_diameter(method=meth, units='microns')


        # Compute RV error #
        rv_std_all = (3 * 10**8) * (center_std) / ((150000)*(diameter))
        if per_10:
            rv_std_avg = (3 * 10**8) * (center_avg_std) / ((150000)*(diameter_avg))
        else:
            rv_std_avg = (3 * 10**8) * (center_avg_std) / ((150000)*(diameter))

        # Convert to m/s #
        center_ms = [(3 * 10**8) * (x) / ((150000)*(diameter)) for x in center]
        center_avg_ms = [(3 * 10**8) * (x) / ((150000)*(diameter)) for x in center_avg]

        # Plot #
        center_line = plt.plot(center_ms, color='g', label='$\sigma_{rv}=%.2f$' % (rv_std_all))
        avg_line = plt.plot(xrange(0, 300, 30), center_avg_ms, color='r', label='$\sigma_{rv}=%.2f$' % (rv_std_avg))

        plt.ylabel('Center drift (m/s)')
        plt.xlabel('Frame number')
        if meth is 'full':
            if cam is 'nf':
                plt.ylim(-5, 5)
            if cam is 'ff':
                plt.ylim(-2, 2)

        plt.legend(loc='lower right')

        # Save #
        if per_10:
            plt.savefig(folder + 'plots_new_avg/rv_error_plots/%s_%s_%s_rv_error.png' % (cam, meth, num_images), bbox_inches='tight')
            print('Saved figure to %splots_new_avg/rv_error_plots/' % str(folder))
        else:
            plt.savefig(folder + 'plots/rv_error_plots/%s_%s_%s_rv_error.png' % (cam, meth, num_images), bbox_inches='tight')
            print('Saved figure to %splots/rv_error_plots/' % str(folder))
        plt.close()

if __name__ == '__main__':
    find_rv_error()
