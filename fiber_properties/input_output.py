"""InputOutput.py was written by Ryan Petersburg for use with fiber
characterization on the EXtreme PREcision Spectrograph
"""
import os
import cPickle as pickle
from collections import Iterable
import matplotlib.pyplot as plt
from astropy.io import fits

def image_list(image_name, ext='.fit', num=10):
    """List of images typically created by FCS."""
    return [image_name + str(i).zfill(3) + ext for i in xrange(num)]

def save_image(input_array, save_file):
    """Saves a np.ndarry as the designated file."""
    save_file = true_file_name(save_file)
    create_directory(save_file)
    # if not (save_file.startswith('C:/') or save_file.startswith('/')
    #         or save_file.startswith('./') or save_file.startswith('../')):
    #     save_file = './' + save_file
    if os.path.split(save_file)[1] in os.listdir(os.path.split(save_file[0])):
    # if save_file.split('/')[-1] in os.listdir('/'.join(save_file.split('/')[:-1])):
        os.remove(save_file)
    if save_file[-3:] in ['tif', 'png', 'pdf', 'eps', 'svg']:
        plt.imsave(save_file, input_array, cmap='gray')
    elif save_file[-3:] == 'fit':
        fits.PrimaryHDU(input_array).writeto(save_file)
    else:
        raise RuntimeError('Please choose either .fit or other standard image extension for file extension')

def save_image_object(image_obj, file_name):
    """Pickle an ImageAnalysis object to file_name."""
    file_name = true_file_name(file_name)
    if file_name[-2:] != '.p' and file_name[-4:] != '.pkl':
        raise RuntimeError('Please use .p or .pkl for file extension')
    create_directory(file_name)
    with open(file_name, 'wb') as output_file:
        pickle.dump(image_obj, output_file, -1)

def load_image_object(object_file, image_file=None):
    """Load a pickled ImageAnalysis object."""
    if object_file[-2:] != '.p' and object_file[-4:] != '.pkl':
        raise RuntimeError('Please use .p or .pkl for file extension')
    with open(object_file, 'rb') as input_file:
        image_obj = pickle.load(input_file)
    if image_file is not None:
        image_obj.set_image_file(image_file)
    return image_obj

def create_directory(file_name):
    """Recursively creates directories if they don't exist."""
    file_name = true_file_name(file_name)
    file_list = file_name.split(os.sep)
    if '.' in file_list[-1]:
        file_list = file_list[:-1]
    directory = os.sep.join(file_list)
    if not os.path.exists(directory):
        print 'Making directory', directory
        os.makedirs(directory)

    # if not (file_name.startswith('C:/') or file_name.startswith('/')
    #         or file_name.startswith('./') or file_name.startswith('../')):
    #     file_name = './' + file_name
    # file_list = file_name.split('/')

    # if '.' in file_list[-1]:
    #     length = len(file_list) - 2
    # else:
    #     length = len(file_list) - 1

    # for i in xrange(length):
    #     if file_list[i+1] not in os.listdir('/'.join(file_list[:i+1])+'/'):
    #         print 'Making directory', '/'.join(file_list[:i+2])
    #         os.mkdir('/'.join(file_list[:i+2]))

def save_data(image_obj, file_name):
    """Save object data as a dictionary in a text file."""
    file_name = true_file_name(file_name)
    data = to_dict(image_obj)
    if file_name[-3:] == 'txt':
        create_directory(file_name)
        with open(file_name, 'w') as save_file:
            save_file.write(str(data))
    else:
        raise RuntimeError('Please use .txt for file extension')

def to_dict(obj):
    """Recursively convert a Python object graph to a dictionary"""
    if isinstance(obj, basestring):
        return obj
    elif isinstance(obj, dict):
        return dict((key, to_dict(val)) for key, val in obj.items())
    elif isinstance(obj, Iterable):
        return [to_dict(val) for val in obj]
    elif hasattr(obj, '__dict__'):
        return to_dict(vars(obj))
    return obj

def true_file_name(file_name):
    return os.path.realpath(file_name)
