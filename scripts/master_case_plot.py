from fiber_properties import FiberImage
import numpy as np
import matplotlib.pyplot as plt

FOLDER = '/Users/Dominic/Box Sync/Fiber_Characterization/Image Analysis/data/modal_noise/rv_error/'

CASES = ['coupled_agitation/', 'LED/', 'slow_agitation/', 'coupled_ag_new/']
# List of 4 lists - one for each case - [[nf, ff], [nf, ff]...]
ANGLE = [[np.deg2rad(90-59), np.deg2rad(90+59)], [np.deg2rad(90-58.5), np.deg2rad(90+58.5)], [np.deg2rad(90-58), np.deg2rad(90+58)], [np.deg2rad(90-52.5), np.deg2rad(90+52.5)]]

METHOD = 'full'
CAMERAS = ['nf']

def case_plot(cases=CASES, folder=FOLDER, camera=CAMERAS, meth=METHOD, angle=ANGLE):

    for cam in camera:
        # Plot figure #
        fig1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, figsize=[10,10])
        fig2, (axis1, axis2, axis3) = plt.subplots(3, 1, sharex=True, figsize=[10,8])
        # Grid spacing #
        fig1.subplots_adjust(hspace=0)
        fig2.subplots_adjust(hspace=0)
        # Labels #
        ax1.set_ylabel('Center drift (m/s)')
        ax3.set_ylabel('Center drift (m/s)')
        ax3.set_xlabel('Frame number')
        ax4.set_xlabel('Frame number')
        axis1.set_ylabel('Center drift (m/s)')
        axis2.set_ylabel('Center drift (m/s)')
        axis3.set_ylabel('Center drift (m/s)')
        axis3.set_xlabel('Frame Number')
        print('Camera: %s' % cam)
        for case in cases:
            print('### Case: %s ###' % case)
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
            rv_std_all = (3 * 10**8) * (center_std) / ((50000)*(diameter))
            rv_std_10avg = (3 * 10**8) * (center_10avg_std) / ((50000)*(diameter))
            rv_std_30avg = (3 * 10**8) * (center_30avg_std) / ((50000)*(diameter))

            # Convert to m/s #
            center_ms = [(3 * 10**8) * (x) / ((50000)*(diameter)) for x in center]
            center_10avg_ms = [(3 * 10**8) * (x) / ((50000)*(diameter)) for x in center_avg_10]
            center_30avg_ms = [(3 * 10**8) * (x) / ((50000)*(diameter)) for x in center_avg_30]

            # Plot #
            linecol1 = 'orange'
            linecol2 = 'red'
            linecol3 = 'black'
            ls1 = '-'
            ls2 = '--'
            ls3 = '-'
            bbox = [0, -0.03, 1, 1]
            fntsize = 11
            if case == 'coupled_agitation/':
                print('plotting case: coupled_agitation')
                ax1.text(0.5, 0.95, 'Coupled Agitation', horizontalalignment='center', transform=ax1.transAxes, fontweight='bold')
                ax1.set_ylabel('Center drift (m/s)')
                center_line, = ax1.plot(center_ms, linecol1, ls=ls1, label='1s: $\sigma_{rv}=%.2f$' % (rv_std_all))
                avg10_line, = ax1.plot(xrange(5, 300, 10), center_10avg_ms, linecol2, ls=ls2, linewidth=1.5, label='10s: $\sigma_{rv}=%.2f$' % (rv_std_10avg))
                avg30_line, = ax1.plot(xrange(15, 300, 30), center_30avg_ms, linecol3, ls=ls3, label='30s: $\sigma_{rv}=%.2f$' % (rv_std_30avg))
                ax1.legend(loc='lower right', fontsize=fntsize, fancybox=True, framealpha=0.3)
            if case == 'LED/':
                print('plotting case: LED')
                ax2.text(0.5, 0.95, 'LED', horizontalalignment='center', transform=ax2.transAxes, fontweight='bold')
                ax2.set_ylabel('Center drift (m/s)')
                center_line, = ax2.plot(center_ms, color=linecol1, ls=ls1, label='1s: $\sigma_{rv}=%.2f$' % (rv_std_all))
                avg10_line, = ax2.plot(xrange(5, 300, 10), center_10avg_ms, linecol2, ls=ls2, linewidth=1.5, label='10s: $\sigma_{rv}=%.2f$' % (rv_std_10avg))
                avg30_line, = ax2.plot(xrange(15, 300, 30), center_30avg_ms, linecol3, ls=ls3, label='30s: $\sigma_{rv}=%.2f$' % (rv_std_30avg))
                ax2.legend(loc='lower right', fontsize=fntsize, fancybox=True, framealpha=0.3)

                axis3.text(0, 0.9, 'LED', horizontalalignment='left', transform=axis3.transAxes, fontweight='bold')
                center_line, = axis3.plot(center_ms, color=linecol1, ls=ls1, label='1s: $\sigma_{rv}=%.2f$' % (rv_std_all))
                avg10_line, = axis3.plot(xrange(5, 300, 10), center_10avg_ms, linecol2, ls=ls2, linewidth=1.5, label='10s: $\sigma_{rv}=%.2f$' % (rv_std_10avg))
                avg30_line, = axis3.plot(xrange(15, 300, 30), center_30avg_ms, linecol3, ls=ls3, label='30s: $\sigma_{rv}=%.2f$' % (rv_std_30avg))
                axis3.legend(bbox_to_anchor=bbox, ncol=3, mode='tight', loc='lower center', fontsize=fntsize, fancybox=True, framealpha=0)
            if case == 'coupled_ag_new/':
                print('plotting case: coupled_ag_new')
                ax3.text(0.5, 0.95, 'New Coupled Agitation', horizontalalignment='center', transform=ax3.transAxes, fontweight='bold')
                ax3.set_ylabel('Center drift (m/s)')
                center_line, = ax3.plot(center_ms, color=linecol1, ls=ls1, label='1s: $\sigma_{rv}=%.2f$' % (rv_std_all))
                avg10_line, = ax3.plot(xrange(5, 300, 10), center_10avg_ms, linecol2, ls=ls2, linewidth=1.5, label='10s: $\sigma_{rv}=%.2f$' % (rv_std_10avg))
                avg30_line, = ax3.plot(xrange(15, 300, 30), center_30avg_ms, linecol3, ls=ls3, label='30s: $\sigma_{rv}=%.2f$' % (rv_std_30avg))
                ax3.legend(loc='lower right', fontsize=fntsize, fancybox=True, framealpha=0.3)

                axis2.text(0, 0.9, 'Coupled Agitation', horizontalalignment='left', transform=axis2.transAxes, fontweight='bold')
                center_line, = axis2.plot(center_ms, color=linecol1, ls=ls1, label='1s: $\sigma_{rv}=%.2f$' % (rv_std_all))
                avg10_line, = axis2.plot(xrange(5, 300, 10), center_10avg_ms, linecol2, ls=ls2, linewidth=1.5, label='10s: $\sigma_{rv}=%.2f$' % (rv_std_10avg))
                avg30_line, = axis2.plot(xrange(15, 300, 30), center_30avg_ms, linecol3, ls=ls3, label='30s: $\sigma_{rv}=%.2f$' % (rv_std_30avg))
                axis2.legend(bbox_to_anchor=bbox, ncol=3, mode='tight', loc='lower center', fontsize=fntsize, framealpha=0)
            if case == 'slow_agitation/':
                print('plotting case: slow_agitation')
                ax4.text(0.5, 0.95, 'Slow Agitation', horizontalalignment='center', transform=ax4.transAxes, fontweight='bold')
                ax4.set_ylabel('Center drift (m/s)')
                center_line, = ax4.plot(center_ms, color=linecol1, ls=ls1, label='1s: $\sigma_{rv}=%.2f$' % (rv_std_all))
                avg10_line, = ax4.plot(xrange(5, 300, 10), center_10avg_ms, linecol2, ls=ls2, linewidth=1.5, label='10s: $\sigma_{rv}=%.2f$' % (rv_std_10avg))
                avg30_line, = ax4.plot(xrange(15, 300, 30), center_30avg_ms, linecol3, ls=ls3, label='30s: $\sigma_{rv}=%.2f$' % (rv_std_30avg))
                ax4.legend(loc='lower right', fontsize=fntsize, fancybox=True, framealpha=0.3)

                axis1.text(0, 0.9, 'Slow Agitation', horizontalalignment='left', transform=axis1.transAxes, fontweight='bold')
                center_line, = axis1.plot(center_ms, color=linecol1, ls=ls1, label='1s: $\sigma_{rv}=%.2f$' % (rv_std_all))
                avg10_line, = axis1.plot(xrange(5, 300, 10), center_10avg_ms, linecol2, ls=ls2, linewidth=1.5, label='10s: $\sigma_{rv}=%.2f$' % (rv_std_10avg))
                avg30_line, = axis1.plot(xrange(15, 300, 30), center_30avg_ms, linecol3, ls=ls3, label='30s: $\sigma_{rv}=%.2f$' % (rv_std_30avg))
                axis1.legend(bbox_to_anchor=bbox, ncol=3, mode='tight', loc='lower center', fontsize=fntsize, fancybox=True, framealpha=0)

        # Finish plotting #
        fig1.tight_layout()
        if cam is 'nf':

            axis1.set_ylim(-40, 40)
            axis1.set_yticks(np.arange(-30, 31, 15))
            axis2.set_ylim(-40, 40)
            axis2.set_yticks(np.arange(-30, 31, 15))
            axis3.set_ylim(-40, 40)
            axis3.set_yticks(np.arange(-30, 31, 15))

        # Save #
        fig1.savefig(folder + 'case_plots/%s_%s_all4_rv_error.png' % (cam, meth), bbox_inches='tight')
        fig2.savefig(folder + 'case_plots/%s_%s_SCL_rv_error.png' % (cam, meth), bbox_inches='tight')

if __name__ == '__main__':
    case_plot()
