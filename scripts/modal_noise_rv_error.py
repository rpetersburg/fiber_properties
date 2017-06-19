from fiber_properties import FiberImage
import matplotlib.pyplot as plt
import os

NEW_DATA = True
NUM_IMAGES = 10
CASE = 1
FOLDER = 'C:/Libraries/Box Sync/ExoLab/Fiber_Characterization/Image Analysis/data/modal_noise/rv_error/'
CAMERAS = ['nf', 'ff']

if CASE == 1:
    FOLDER += 'coupled_agitation/'
if CASE == 2:
    FOLDER += 'LED/'

def main(folder=FOLDER, cameras=CAMERAS, num_images=NUM_IMAGES, new_data=NEW_DATA):
    for cam in cameras:
        centers = []
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
            center = im_obj.get_fiber_center(method='circle', units='microns') - im_obj.get_fiber_centroid(method='full', units='microns')
            print center
            centers.append(center)
            im_obj.save_object(object_file)

        plt.figure()
        plt.subplot(211)
        plt.plot([center.x for center in centers])

        plt.subplot(212)
        plt.plot([center.y for center in centers])

        plt.show()

if __name__ == '__main__':
    main()
