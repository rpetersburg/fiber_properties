from fiber_properties import FiberImage
import numpy as np

# Script to transfer the calculated fiber center and centroid to from one pkl file set to another #

# Data is copied from the FILEPATH to the TARGET files #
# METHOD is the method for which you want to copy to the target files #

FILEPATH = '/Users/Dominic/Box Sync/Fiber_Characterization/Image Analysis/data/modal_noise/rv_error/'
TARGET = '/Users/Dominic/Box Sync/Fiber_Characterization/Image Analysis/data/modal_noise/rv_error/'
METHOD = 'full'

CASE = 1
NUM_IMAGES = 1
CAMERAS = ['ff']

if CASE == 1:
    FILEPATH += 'coupled_agitation/full_pkls/'
    TARGET += 'coupled_agitation/'
if CASE == 2:
    FILEPATH += 'LED/full_pkls/'
    TARGET += 'LED/'
if CASE == 3:
    FILEPATH += 'slow_agitation/full_pkls/'
    TARGET += 'slow_agitation/'


def condenser(filepath=FILEPATH, target=TARGET, meth=METHOD, camera=CAMERAS, num_images=NUM_IMAGES):
    for cam in camera:
        print('Setting fiber center and centroid for %s objects' % cam)
        for i in xrange(0, 300, num_images):
            # For center #
            target_obj = target + cam + '_' + str(i).zfill(3) + '_obj.pkl'
            im_obj = FiberImage(target_obj)
            x_center = getattr(im_obj._center, meth).x
            y_center = getattr(im_obj._center, meth).y

            if type(x_center) is not np.float64:
                if type(y_center) is not np.float64:  # If target_obj does not have data
                    # Copy data from FILEPATH object to TARGET object #
                    copy_obj = filepath + cam + '_' + str(i).zfill(3) + '_obj.pkl'
                    im_obj2 = FiberImage(copy_obj)
                    set_value = getattr(im_obj2._center, meth)
                    setattr(im_obj._center, meth, set_value)
                    im_obj.save_object(target_obj)

            if type(x_center) is np.float64:
                if type(y_center) is np.float64 or float:
                    print('Object already has center data for %s method' % meth)

            # For centroid #
            x_centroid = getattr(im_obj._centroid, meth).x
            y_centroid = getattr(im_obj._centroid, meth).y

            if type(x_centroid) is not np.float64:
                if type(y_centroid) is not np.float64:
                    copy_obj = filepath + cam + '_' + str(i).zfill(3) + '_obj.pkl'
                    im_obj2 = FiberImage(copy_obj)
                    set_value = getattr(im_obj2._centroid, meth)
                    setattr(im_obj._centroid, meth, set_value)
                    im_obj.save_object(target_obj)

            if type(x_centroid) is np.float64:
                if type(y_centroid) is np.float64 or float:
                    print('Object already has centroid data for %s method' % meth)

if __name__ == '__main__':
    condenser()