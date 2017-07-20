from fiber_properties import FiberImage
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

FOLDER = '/Users/Dominic/Box Sync/Fiber_Characterization/Image Analysis/data/modal_noise/rv_error/'

CASES = ['coupled_agitation/', 'LED/', 'slow_agitation/', 'coupled_ag_new/']
# List of 4 lists - one for each case - [[nf, ff], [nf, ff]...]
ANGLE = [[np.deg2rad(90-59), np.deg2rad(90+59)], [np.deg2rad(90-58.5), np.deg2rad(90+58.5)], [np.deg2rad(90-58), np.deg2rad(90+58)], [np.deg2rad(90-52.5), np.deg2rad(90+52.5)]]

METHOD = 'full'
CAMERAS = ['nf', 'ff']

def case_plot(cases=CASES, folder=FOLDER, camera=CAMERAS, meth=METHOD, angle=ANGLE):

    for cam in camera:
        # Plot figure #
        fig1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, figsize=[10,10])
        fig2, (axis1, axis2, axis3) = plt.subplots(3, 1, sharex=True, figsize=[10,10])
        ax1.set_ylabel('Center drift (m/s)')
        ax3.set_ylabel('Center drift (m/s)')
        ax3.set_xlabel('Frame number')
        ax4.set_xlabel('Frame number')
        axis1.set_ylabel('Center drift (m/s)')
        axis2.set_ylabel('Center drift (m/s)')
        axis3.set_ylabel('Center drift (m/s)')
        axis3.set_xlabel('Frame Number')
        for case in cases:
            print('Camera: %s' % cam)
            print('Case: %s' % case)
            print('Method: %s' % meth)
            nf_ang = angle[cases.index(case)][0]
            ff_ang = angle[cases.index(case)][1]
            center_x = []
            center_y = []
            for i in xrange(0, 300, 1):
                im_obj = FiberImage(folder + case + cam + '_' + str(i).zfill(3) + '_obj.pkl')
                c = im_obj.get_fiber_center(method=meth, units='microns') - im_obj.get_fiber_centroid(method=meth, units='microns')
                center_x.append(c.x)
                center_y.append(c.y)
            print('Found fiber centers')

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
                    calc = np.sqrt(c_x**2 + c_y**2) * np.cos(np.arctan(c_y/c_x) + (nf_ang) + (np.pi/2)*(1 - np.sign(c_x)))
                    center.append(calc)
                if cam is 'ff':
                    calc = np.sqrt(c_x**2 + c_y**2) * np.cos(np.arctan(c_y/c_x) + (ff_ang) + (np.pi/2)*(1 - np.sign(c_x)))
                    center.append(calc)

            # Make average lines #
            # 10
            center_avg_10 = []
            itr = range(10)
            for i in xrange(0, 300, 10):
                center_to_avg = []
                for num in itr:
                    center_to_avg.append(center[num])
                center_avg_10.append(np.average(center_to_avg))
                itr = [x + 10 for x in itr]

            # 30 #
            center_avg_30 = []
            itr = range(30)
            for i in xrange(0, 300, 30):
                center_to_avg = []
                for num in itr:
                    center_to_avg.append(center[num])
                center_avg_30.append(np.average(center_to_avg))
                itr = [x + 30 for x in itr]

            # Compute std #
            center_std = np.std(center)
            center_10avg_std = np.std(center_avg_10)
            center_30avg_std = np.std(center_avg_30)

            # Get diameter #
            if cam is 'nf':
                diameter = 100
                diameter_avg = 100
            if cam is 'ff':
                print('Getting ff diameter')
                diameter = im_obj.get_fiber_diameter(method=meth, units='microns')


            # Compute RV error #
            rv_std_all = (3 * 10**8) * (center_std) / ((150000)*(diameter))
            rv_std_10avg = (3 * 10**8) * (center_10avg_std) / ((150000)*(diameter))
            rv_std_30avg = (3 * 10**8) * (center_30avg_std) / ((150000)*(diameter))

            # Convert to m/s #
            center_ms = [(3 * 10**8) * (x) / ((150000)*(diameter)) for x in center]
            center_10avg_ms = [(3 * 10**8) * (x) / ((150000)*(diameter)) for x in center_avg_10]
            center_30avg_ms = [(3 * 10**8) * (x) / ((150000)*(diameter)) for x in center_avg_30]

            # Plot #
            dashes = [25, 5, 25, 5]  # Custom dashes: 10 points on, 5 off, 100 on, 5 off
            if case == 'coupled_agitation/':
                print('plotting case: coupled_agitation')
                ax1.set_title('Coupled Agitation')
                ax1.set_ylabel('Center drift (m/s)')
                center_line, = ax1.plot(center_ms, color='b', label='$\sigma_{rv}=%.2f$' % (rv_std_all))
                avg10_line, = ax1.plot(xrange(5, 300, 10), center_10avg_ms, 'g', label='$\sigma_{rv10}=%.2f$' % (rv_std_10avg), linewidth=2.5)
                avg10_line.set_dashes(dashes)
                avg30_line, = ax1.plot(xrange(15, 300, 30), center_30avg_ms, 'r--', label='$\sigma_{rv30}=%.2f$' % (rv_std_30avg), linewidth=4.0)
                ax1.legend(loc='lower right')
            if case == 'LED/':
                print('plotting case: LED')
                ax2.set_title('LED')
                ax2.set_ylabel('Center drift (m/s)')
                center_line, = ax2.plot(center_ms, color='b', label='$\sigma_{rv}=%.2f$' % (rv_std_all))
                avg10_line, = ax2.plot(xrange(5, 300, 10), center_10avg_ms, 'g', label='$\sigma_{rv10}=%.2f$' % (rv_std_10avg), linewidth=2.5)
                avg10_line.set_dashes(dashes)
                avg30_line, = ax2.plot(xrange(15, 300, 30), center_30avg_ms, 'r--', label='$\sigma_{rv30}=%.2f$' % (rv_std_30avg), linewidth=4.0)
                ax2.legend(loc='lower right')

                axis3.set_title('LED')
                center_line, = axis3.plot(center_ms, color='b', label='$\sigma_{rv}=%.2f$' % (rv_std_all))
                avg10_line, = axis3.plot(xrange(5, 300, 10), center_10avg_ms, 'g', label='$\sigma_{rv10}=%.2f$' % (rv_std_10avg), linewidth=2.5)
                avg10_line.set_dashes(dashes)
                avg30_line, = axis3.plot(xrange(15, 300, 30), center_30avg_ms, 'r--', label='$\sigma_{rv30}=%.2f$' % (rv_std_30avg), linewidth=4.0)
                axis3.legend(loc='lower right')
            if case == 'coupled_ag_new/':
                print('plotting case: coupled_ag_new')
                ax3.set_title('New Coupled Agitation')
                ax3.set_ylabel('Center drift (m/s)')
                center_line, = ax3.plot(center_ms, color='b', label='$\sigma_{rv}=%.2f$' % (rv_std_all))
                avg10_line, = ax3.plot(xrange(5, 300, 10), center_10avg_ms, 'g', label='$\sigma_{rv10}=%.2f$' % (rv_std_10avg), linewidth=2.5)
                avg10_line.set_dashes(dashes)
                avg30_line, = ax3.plot(xrange(15, 300, 30), center_30avg_ms, 'r--', label='$\sigma_{rv30}=%.2f$' % (rv_std_30avg), linewidth=4.0)
                ax3.legend(loc='lower right')

                axis2.set_title('Coupled Agitation')
                center_line, = axis2.plot(center_ms, color='b', label='$\sigma_{rv}=%.2f$' % (rv_std_all))
                avg10_line, = axis2.plot(xrange(5, 300, 10), center_10avg_ms, 'g', label='$\sigma_{rv10}=%.2f$' % (rv_std_10avg), linewidth=2.5)
                avg10_line.set_dashes(dashes)
                avg30_line, = axis2.plot(xrange(15, 300, 30), center_30avg_ms, 'r--', label='$\sigma_{rv30}=%.2f$' % (rv_std_30avg), linewidth=4.0)
                axis2.legend(loc='lower right')
            if case == 'slow_agitation/':
                print('plotting case: slow_agitation')
                ax4.set_title('Slow Agitation')
                ax4.set_ylabel('Center drift (m/s)')
                center_line, = ax4.plot(center_ms, color='b', label='$\sigma_{rv}=%.2f$' % (rv_std_all))
                avg10_line, = ax4.plot(xrange(5, 300, 10), center_10avg_ms, 'g', label='$\sigma_{rv10}=%.2f$' % (rv_std_10avg), linewidth=2.5)
                avg10_line.set_dashes(dashes)
                avg30_line, = ax4.plot(xrange(15, 300, 30), center_30avg_ms, 'r--', label='$\sigma_{rv30}=%.2f$' % (rv_std_30avg), linewidth=4.0)
                ax4.legend(loc='lower right')

                axis1.set_title('Slow Agitation')
                center_line, = axis1.plot(center_ms, color='b', label='$\sigma_{rv}=%.2f$' % (rv_std_all))
                avg10_line, = axis1.plot(xrange(5, 300, 10), center_10avg_ms, 'g', label='$\sigma_{rv10}=%.2f$' % (rv_std_10avg), linewidth=2.5)
                avg10_line.set_dashes(dashes)
                avg30_line, = axis1.plot(xrange(15, 300, 30), center_30avg_ms, 'r--', label='$\sigma_{rv30}=%.2f$' % (rv_std_30avg), linewidth=4.0)
                axis1.legend(loc='lower right')

        # Finish plotting #
        fig1.tight_layout()
        if cam is 'nf':
            ax1.set_ylim(-5, 5)
            ax2.set_ylim(-5, 5)
            ax3.set_ylim(-5, 5)
            ax4.set_ylim(-20, 20)

            axis1.set_ylim(-20, 20)
            axis2.set_ylim(-5, 5)
            axis3.set_ylim(-5, 5)
        if cam is 'ff':
            ax1.set_ylim(-2, 2)
            ax2.set_ylim(-2, 2)
            ax3.set_ylim(-2, 2)
            ax4.set_ylim(-15, 15)

            axis1.set_ylim(-15, 15)
            axis2.set_ylim(-2, 2)
            axis3.set_ylim(-2, 2)

        # Save #
        fig1.savefig(folder + 'case_plots/%s_%s_all4_rv_error.png' % (cam, meth), bbox_inches='tight')
        fig2.savefig(folder + 'case_plots/%s_%s_SCL_rv_error.png' % (cam, meth), bbox_inches='tight')

if __name__ == '__main__':
    case_plot()
