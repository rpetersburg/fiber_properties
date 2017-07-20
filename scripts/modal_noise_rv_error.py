from fiber_properties import FiberImage
import matplotlib.pyplot as plt
import os
import numpy as np
from multiprocessing import Pool

PLOTTING = True
MULTIPROCESS = True
PROCESSES = 6

PLOT_PER_NUM = True  # New method of plotting per 10 average over NUM_IMAGES = 1
NUMBER = 10

PLOT_FIBER_CENTROID = True
NEW_DATA = False
NUM_IMAGES = 1
CASE = 1
CAMERAS = ['nf', 'ff']

FOLDER = '/Users/Dominic/Box Sync/Fiber_Characterization/Image Analysis/data/modal_noise/rv_error/'
METHOD = 'full'
CENTER_RANGE = None  # for circle method. Default None

if CASE == 1:
    FOLDER += 'coupled_agitation/'
if CASE == 2:
    FOLDER += 'LED/'
if CASE == 3:
    FOLDER += 'slow_agitation/'
if CASE == 4:
    FOLDER += 'coupled_ag_new/'


def multi(processes=PROCESSES, camera=CAMERAS, num_images=NUM_IMAGES):
    p = Pool(processes)
    for cam in camera:
        objects = []
        for i in xrange(0, 300, num_images):
            if num_images == 1:
                object_file = cam + '_' + str(i).zfill(3) + '_obj.pkl'
            else:
                object_file = cam + '_' + str(i).zfill(3) + '-' + str(i+num_images-1).zfill(3) + '_obj.pkl'
            objects.append(object_file)
        p.map(fiber_center_multi, objects)


def fiber_center_multi(object_file, folder=FOLDER, cameras=CAMERAS, num_images=NUM_IMAGES, new_data=NEW_DATA, meth=METHOD, plot_fiber_centroid=PLOT_FIBER_CENTROID, cent_range=CENTER_RANGE):
    cam = object_file[:2]
    ject = object_file[:6]
    if object_file not in os.listdir(folder) or new_data:
        images = [folder + ject + '.fit']
        ambient = [folder + 'ambient/' + cam + '_' + str(j).zfill(3) + '.fit' for j in xrange(10)]
        dark = [folder + 'dark/' + cam + '_' + str(j).zfill(3) + '.fit' for j in xrange(10)]
        im_obj = FiberImage(images, ambient=ambient, dark=dark, camera=cam)
        im_obj.save_object(folder + object_file)

    object_file = folder + object_file
    im_obj = FiberImage(object_file)

    print('Getting fiber center for image %s...' % object_file[-14:])
    fiber_centroid = im_obj.get_fiber_center(method=meth, units='microns', center_range=cent_range) - im_obj.get_fiber_centroid(method=meth, units='microns', center_range=cent_range)
    im_obj.save_object(object_file)


def fiber_center(folder=FOLDER, cameras=CAMERAS, num_images=NUM_IMAGES, new_data=NEW_DATA, meth=METHOD, plot_fiber_centroid=PLOT_FIBER_CENTROID, per_num=PLOT_PER_NUM, number=NUMBER, cent_range=CENTER_RANGE):
    for cam in cameras:
        print('Saving to Folder: %s' % folder)
        print('Plotting fiber centroid: %s' % str(plot_fiber_centroid))
        print('Camera: %s' % str(cam))
        print('Method: %s' % str(meth))
        if plot_fiber_centroid:
            all_x = []
            all_y = []
        else:
            center_x = []
            center_y = []
            centroid_x = []
            centroid_y = []
        for i in xrange(0, 300, num_images):
            if num_images == 1:
                object_file = cam + '_' + str(i).zfill(3) + '_obj.pkl'
            else:
                object_file = cam + '_' + str(i).zfill(3) + '-' + str(i+num_images-1).zfill(3) + '_obj.pkl'

            if object_file not in os.listdir(folder) or new_data:
                images = [folder + cam + '_' + str(j).zfill(3) + '.fit' for j in xrange(i, i+num_images)]
                ambient = [folder + 'ambient/' + cam + '_' + str(j).zfill(3) + '.fit' for j in xrange(10)]
                dark = [folder + 'dark/' + cam + '_' + str(j).zfill(3) + '.fit' for j in xrange(10)]
                im_obj = FiberImage(images, ambient=ambient, dark=dark, camera=cam)
                im_obj.save_object(folder + object_file)

            object_file = folder + object_file
            im_obj = FiberImage(object_file)

            if plot_fiber_centroid:
                print('Getting fiber center for image %s ...' % str(i))
                fiber_centroid = im_obj.get_fiber_center(method=meth, units='microns', center_range=cent_range) - im_obj.get_fiber_centroid(method=meth, units='microns', center_range=cent_range)
                all_x.append(float(fiber_centroid.x))
                all_y.append(float(fiber_centroid.y))
                print(fiber_centroid)
            else:
                print('Getting center and centroid for image %s ...' % str(i))
                center = im_obj.get_fiber_center(method=meth, units='microns')
                centroid = im_obj.get_fiber_centroid(method=meth, units='microns')
                center_x.append(float(center.x))
                center_y.append(float(center.y))
                centroid_x.append(float(centroid.x))
                centroid_y.append(float(centroid.y))
                print('center: %s' % center)
                print('centroid: %s' % centroid)

            im_obj.save_object(object_file)

        ##########

        # Get drift and std #
        if plot_fiber_centroid:
            x_median = np.median(all_x)
            y_median = np.median(all_y)
            drift_x = [x - x_median for x in all_x]
            drift_y = [y - y_median for y in all_y]

            std_x = np.std(drift_x)
            std_y = np.std(drift_y)
            sig_xy = np.sqrt((std_x**2)+(std_y**2))

        else:
            center_x_median = np.median(center_x)
            center_y_median = np.median(center_y)
            centroid_x_median = np.median(centroid_x)
            centroid_y_median = np.median(centroid_y)
            center_drift_x = [x - center_x_median for x in center_x]
            center_drift_y = [y - center_y_median for y in center_y]
            centroid_drift_x = [x - centroid_x_median for x in centroid_x]
            centroid_drift_y = [y - centroid_y_median for y in centroid_y]

            center_std_x = np.std(center_drift_x)
            center_std_y = np.std(center_drift_y)
            center_sig_xy = np.sqrt((center_std_x**2)+(center_std_y**2))
            centroid_std_x = np.std(centroid_drift_x)
            centroid_std_y = np.std(centroid_drift_y)
            centroid_sig_xy = np.sqrt((centroid_std_x**2)+(centroid_std_y**2))

        ##########

        # Get average per num frames #
        if per_num:
            if plot_fiber_centroid:
                avg_xlist = []
                avg_ylist = []
            else:
                center_xavg = []
                center_yavg = []
                centroid_xavg = []
                centroid_yavg = []
            for i in xrange(0, 300, number):
                avg_file = cam + '_' + str(i).zfill(3) + '-' + str(i+num-1).zfill(3) + '_obj.pkl'

                if avg_file not in os.listdir(folder) or new_data:
                    images = [folder + cam + '_' + str(j).zfill(3) + '.fit' for j in xrange(i, i+number)]
                    ambient = [folder + 'ambient/' + cam + '_' + str(j).zfill(3) + '.fit' for j in xrange(10)]
                    dark = [folder + 'dark/' + cam + '_' + str(j).zfill(3) + '.fit' for j in xrange(10)]
                    avg_obj = FiberImage(images, ambient=ambient, dark=dark, camera=cam)
                    avg_obj.save_object(folder + avg_file)

                avg_file = folder + avg_file
                avg_obj = FiberImage(avg_file)

                if plot_fiber_centroid:
                    print('Getting fiber center for image %s ...' % str(i))
                    fiber_centroid = avg_obj.get_fiber_center(method=meth, units='microns', center_range=cent_range) - avg_obj.get_fiber_centroid(method=meth, units='microns', center_range=cent_range)
                    avg_xlist.append(float(fiber_centroid.x))
                    avg_ylist.append(float(fiber_centroid.y))
                    print(fiber_centroid)
                else:
                    print('Getting center and centroid for image %s ...' % str(i))
                    center = avg_obj.get_fiber_center(method=meth, units='microns')
                    centroid = avg_obj.get_fiber_centroid(method=meth, units='microns')
                    center_xavg.append(float(center.x))
                    center_yavg.append(float(center.y))
                    centroid_xavg.append(float(centroid.x))
                    centroid_yavg.append(float(centroid.y))
                    print('center: %s' % center)
                    print('centroid: %s' % centroid)

        else:
            if plot_fiber_centroid:
                avg_xlist = []
                avg_ylist = []
                itr = range(number)
                for i in xrange(0, 300, number):
                    x_to_avg = []
                    y_to_avg = []
                    for num in itr:
                        x_to_avg.append(drift_x[num])
                        y_to_avg.append(drift_y[num])
                    avg_xlist.append(np.average(x_to_avg))
                    avg_ylist.append(np.average(y_to_avg))
                    itr = [x + number for x in itr]

                std_x_avg = np.std(avg_xlist)
                std_y_avg = np.std(avg_ylist)
                sig_xy_avg = np.sqrt((std_x_avg**2)+(std_y_avg**2))

            else:
                center_xavg = []
                center_yavg = []
                centroid_xavg = []
                centroid_yavg = []
                itr = range(number)
                for i in xrange(0, 300, number):
                    cenx_to_avg = []
                    ceny_to_avg = []
                    centx_to_avg = []
                    centy_to_avg = []
                    for num in itr:
                        cenx_to_avg.append(center_drift_x[num])
                        ceny_to_avg.append(center_drift_y[num])
                        centx_to_avg.append(centroid_drift_x[num])
                        centy_to_avg.append(centroid_drift_y[num])
                    center_xavg.append(np.average(cenx_to_avg))
                    center_yavg.append(np.average(ceny_to_avg))
                    centroid_xavg.append(np.average(centx_to_avg))
                    centroid_yavg.append(np.average(centy_to_avg))
                    itr = [x + number for x in itr]

                center_std_x_avg = np.std(center_xavg)
                center_std_y_avg = np.std(center_yavg)
                center_sig_xy_avg = np.sqrt((center_std_x_avg**2)+(center_std_y_avg**2))
                centroid_std_x_avg = np.std(centroid_xavg)
                centroid_std_y_avg = np.std(centroid_yavg)
                centroid_sig_xy_avg = np.sqrt((centroid_std_x_avg**2)+(centroid_std_y_avg**2))

        ##########

        # Get drift and std of avg #
        if per_num:
            if plot_fiber_centroid:
                xavg_median = np.median(avg_xlist)
                yavg_median = np.median(avg_ylist)
                avg_xlist = [x - x_median for x in avg_xlist]
                avg_ylist = [y - y_median for y in avg_ylist]

                std_x_avg = np.std(avg_xlist)
                std_y_avg = np.std(avg_ylist)
                sig_xy_avg = np.sqrt((std_x_avg**2)+(std_y_avg**2))

            else:
                center_xavg_median = np.median(center_xavg)
                center_yavg_median = np.median(center_yavg)
                centroid_xavg_median = np.median(centroid_xavg)
                centroid_yavg_median = np.median(centroid_yavg)
                center_xavg = [x - center_xavg_median for x in center_xavg]
                center_yavg = [y - center_yavg_median for y in center_yavg]
                centroid_xavg = [x - centroid_xavg_median for x in centroid_xavg]
                centroid_yavg = [y - centroid_yavg_median for y in centroid_yavg]

                center_std_x_avg = np.std(center_xavg)
                center_std_y_avg = np.std(center_yavg)
                center_sig_xy_avg = np.sqrt((center_std_x_avg**2)+(center_std_y_avg**2))
                centroid_std_x_avg = np.std(centroid_xavg)
                centroid_std_y_avg = np.std(centroid_yavg)
                centroid_sig_xy_avg = np.sqrt((centroid_std_x_avg**2)+(centroid_std_y_avg**2))

        ##########

        if PLOTTING:
            # Plot #
            if plot_fiber_centroid:
                fig1, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
                drift_xline, = ax1.plot(drift_x, 'b')
                drift_yline, = ax2.plot(drift_y, 'b')

                if num_images == 1 or per_num:
                    avg_xline, = ax1.plot(xrange(number/2, 300, number), avg_xlist, 'r')
                    avg_yline, = ax2.plot(xrange(number/2, 300, number), avg_ylist, 'r')

                ax1.set_ylabel('X drift ($\mu m$)')
                ax2.set_ylabel('Y drift ($\mu m$)')
                ax2.set_xlabel('Frame number')

            else:
                fig1, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
                fig2, (ax3, ax4) = plt.subplots(2, 1, sharex=True)
                center_xline, = ax1.plot(center_drift_x, 'b')
                center_yline, = ax2.plot(center_drift_y, 'b')
                centroid_xline, = ax3.plot(centroid_drift_x, 'g')
                centroid_yline, = ax4.plot(centroid_drift_y, 'g')

                if num_images == 1 or per_num:
                    center_avg_xline, = ax1.plot(xrange(number/2, 300, number), center_xavg, 'r')
                    center_avg_yline, = ax2.plot(xrange(number/2, 300, number), center_yavg, 'r')
                    centroid_avg_xline, = ax3.plot(xrange(number/2, 300, number), centroid_xavg, 'c')
                    centroid_avg_yline, = ax4.plot(xrange(number/2, 300, number), centroid_yavg, 'c')

                ax1.set_ylabel('Center X\ndrift ($\mu m$)')
                ax2.set_ylabel('Center Y\ndrift ($\mu m$)')
                ax2.set_xlabel('Frame number')
                ax3.set_ylabel('Centroid X\ndrift ($\mu m$)')
                ax4.set_ylabel('Centroid Y\ndrift ($\mu m$)')
                ax4.set_xlabel('Frame number')

            if plot_fiber_centroid:
                if num_images == 1 or per_num:
                    fig1.legend((drift_xline, avg_xline), ('$\sigma_x=%.2f$\n$\sigma_y=%.2f$\n$\sigma_t=%.2f$' % (std_x, std_y, sig_xy), '$\sigma_x=%.2f$\n$\sigma_y=%.2f$\n$\sigma_t=%.2f$' % (std_x_avg, std_y_avg, sig_xy_avg)), loc='lower center')
                    fig1.subplots_adjust(bottom=0.4)
                else:
                    fig1.legend(drift_xline, '$\sigma_x=%.2f$\n$\sigma_y=%.2f$\n$\sigma_t=%.2f$' % (std_x, std_y, sig_xy), loc='lower center')
                    fig1.subplots_adjust(bottom=0.3)
            else:
                if num_images == 1 or per_num:
                    fig1.legend((center_xline, center_avg_xline), ('$\sigma_x=%.2f$\n$\sigma_y=%.2f$\n$\sigma_t=%.2f$' % (center_std_x, center_std_y, center_sig_xy), '$\sigma_x=%.2f$\n$\sigma_y=%.2f$\n$\sigma_t=%.2f$' % (center_std_x_avg, center_std_y_avg, center_sig_xy_avg)), loc='lower center')
                    fig2.legend((centroid_xline, centroid_avg_xline), ('$\sigma_x=%.2f$\n$\sigma_y=%.2f$\n$\sigma_t=%.2f$' % (centroid_std_x, centroid_std_y, centroid_sig_xy), '$\sigma_x=%.2f$\n$\sigma_y=%.2f$\n$\sigma_t=%.2f$' % (centroid_std_x_avg, centroid_std_y_avg, centroid_sig_xy_avg)), loc='lower center')
                    fig1.subplots_adjust(bottom=0.4)
                    fig2.subplots_adjust(bottom=0.4)
                else:
                    fig1.legend(center_xline, '$\sigma_x=%.2f$\n$\sigma_y=%.2f$\n$\sigma_t=%.2f$' % (center_std_x, center_std_y, center_sig_xy), loc='lower center')
                    fig2.legend(centroid_xline, '$\sigma_x=%.2f$\n$\sigma_y=%.2f$\n$\sigma_t=%.2f$' % (centroid_std_x, centroid_std_y, centroid_std_xy), loc='lower center')
                    fig1.subplots_adjust(bottom=0.3)
                    fig2.subplots_adjust(bottom=0.3)

            # Save #
            if plot_fiber_centroid:
                if per_num:
                    fig1.savefig(folder + 'plots_new_avg/' + '%s_%s_%s_per%savg_img_plot.png' % (cam, num_images, meth, number), bbox_inches='tight')
                    print('Saved figure: %splots_new_avg/%s_%s_%s_per%savg_img_plot.png' % (folder, cam, num_images, meth, number))
                else:
                    fig1.savefig(folder + 'plots/' + '%s_%s_%s_statavg%s_img_plot.png' % (cam, num_images, meth, number), bbox_inches='tight')
                    print('Saved figure: %splots/%s_%s_%s_statav%s_img_plot.png' % (folder, cam, num_images, meth, number))
            else:
                if per_num:
                    fig1.savefig(folder + 'plots_new_avg/' + 'center+centroid/center_%s_%s_%s_per%savg_img_plot.png' % (cam, num_images, meth, number), bbox_inches='tight')
                    fig2.savefig(folder + 'plots_new_avg/' + 'center+centroid/centroid_%s_%s_%s_per%savg_img_plot.png' % (cam, num_images, meth, number), bbox_inches='tight')
                    print('Saved figures to: %splots_new_avg/center+centroid/' % folder)
                else:
                    fig1.savefig(folder + 'plots/' + 'center+centroid/center_%s_%s_%s_statavg%s_img_plot.png' % (cam, num_images, meth, number), bbox_inches='tight')
                    fig2.savefig(folder + 'plots/' + 'center+centroid/centroid_%s_%s_%s_statavg%s_img_plot.png' % (cam, num_images, meth, number), bbox_inches='tight')
                    print('Saved figures to: %splots/center+centroid/' % folder)


if __name__ == '__main__':
    if MULTIPROCESS:
        multi()
        if PLOTTING:
            fiber_center()
    else:
        fiber_center()
