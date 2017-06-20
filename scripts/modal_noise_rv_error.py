from fiber_properties import FiberImage
import matplotlib.pyplot as plt
import os
import numpy as np

NEW_DATA = True
NUM_IMAGES = 10
CASE = 1
FOLDER = '/Users/Dominic/Box Sync/Fiber_Characterization/Image Analysis/data/modal_noise/rv_error/'
CAMERAS = ['nf', 'ff']

if CASE == 1:
    FOLDER += 'coupled_agitation/'
if CASE == 2:
    FOLDER += 'LED/'

def main(folder=FOLDER, cameras=CAMERAS, num_images=NUM_IMAGES, new_data=NEW_DATA):
    for cam in cameras:
        print(folder)
        print('camera: %s' % str(cam))
        centers = []
        all_x = []
        all_y = []
        for i in xrange(0, 300-num_images, num_images):
            object_file = cam + '_' + str(i).zfill(3) + '-' + str(i+num_images).zfill(3) + '_obj.pkl'

            if object_file not in os.listdir(folder) or new_data:
                images = [folder + cam + '_' + str(j).zfill(3) + '.fit' for j in xrange(i, i+num_images)]
                ambient = [folder + 'ambient/' + cam + '_' + str(j).zfill(3) + '.fit' for j in xrange(min(10, num_images))]
                dark = [folder + 'dark/' + cam + '_' + str(j).zfill(3) + '.fit' for j in xrange(min(10, num_images))]
                im_obj = FiberImage(images, ambient=ambient, dark=dark, camera=cam)
                im_obj.save_object(folder + object_file)

            object_file = folder + object_file
            im_obj = FiberImage(object_file)
            print('Getting center for image set %s ...' % str(i))
            center = im_obj.get_fiber_center(method='edge', units='microns') - im_obj.get_fiber_centroid(method='edge', units='microns')
            print(center)
            all_x.append(float(center.x))
            all_y.append(float(center.y))
            im_obj.save_object(object_file)

        drift_x = [x - all_x[0] for x in all_x]
        drift_y = [y - all_y[0] for y in all_y]

        fig1, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        line1 = ax1.plot(drift_x)
        line2, = ax2.plot(drift_y)

        ax1.set_ylabel('X drift ($\mu m$)')
        ax2.set_ylabel('Y drift ($\mu m$)')
        ax2.set_xlabel('Frame number')

        std_x = np.std(all_x)
        std_y = np.std(all_y)
        sig = ((std_x**2)+(std_y**2))**0.5

        fig1.legend(line1, ['$\sigma_x=%.2f$\n$\sigma_y=%.2f$\n$\sigma_t=%.2f$' % (std_x, std_y, sig)], loc='lower center')
        plt.subplots_adjust(bottom=0.3)

        fig1.savefig(folder + 'plots/%s_%s_img_plot.png' % (cam, num_images), bbox_inches='tight')
        print('Saved figure to %splots' % str(folder))

if __name__ == '__main__':
    main()
