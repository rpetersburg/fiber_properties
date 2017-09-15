"""InputOutput.py was written by Ryan Petersburg for use with fiber
characterization on the EXtreme PREcision Spectrograph
"""
import os
# import cPickle as pickle
import pickle
from collections import Iterable
from numbers import Number
import matplotlib.pyplot as plt
from astropy.io import fits
from PIL import Image
import numpy as np

def image_list(image_name, ext='.fit', num=10, start=0):
    """List of images typically created by FCS."""
    return [image_name + str(i).zfill(3) + ext for i in xrange(start, start+num, 1)]

def save_image(input_array, image_file):
    """Saves a np.ndarry as the designated file."""
    image_file = true_path(image_file)
    create_directory(image_file)
    if os.path.split(image_file)[1] in os.listdir(get_directory(image_file)):
        os.remove(image_file)
    if image_file[-3:] in ['tif', 'png', 'pdf', 'eps', 'svg']:
        plt.imsave(image_file, input_array, cmap='gray')
    elif image_file[-3:] == 'fit':
        fits.PrimaryHDU(input_array).writeto(image_file)
    else:
        raise RuntimeError('Please choose either .fit or other standard image'
                           + 'extension for file extension')

def load_image(image_file):
    """Return an image as a numpy array."""
    if image_file[-3:] in ['tif', 'png', 'pdf', 'eps', 'svg']:
        raw_image = Image.open(image_file)
        image = np.array(raw_image).astype('float64')
    elif image_file[-3:] in ['fit']:
        raw_image = fits.open(image_file, ignore_missing_end=True)[0]
        image = raw_image.data.astype('float64')
    return image

def save_header(header, image_file):
    """Save a dictionary to a FITS header."""
    image_header = fits.open(image_file)[0].header
    with fits.open(image_file, 'update') as f:
        for hdu in f:
            image_header = hdu.header
            for key in header:
                if isinstance(header[key], Number):
                    image_header[key] = header[key]
                else:
                    image_header[key] = str(header[key])

def load_header(image_file):
    """Return an image header as a python dictionary."""
    if image_file[-3:] == 'fit':
        raw_image = fits.open(image_file, ignore_missing_end=True)[0]
        header = dict(raw_image.header)
    elif image_file[-3:] == 'tif':
        raw_image = Image.open(image_file)
        # Complicated way to get the header from a TIFF image as a dictionary
        header = dict([i.split('=') for i in raw_image.tag[270][0].split('\r\n')][:-1])
        header['BITPIX'] = int(raw_image.tag[258][0])
    return header

def save_image_object(image_obj, file_name):
    """Pickle an ImageAnalysis object to file_name."""
    file_name = true_path(file_name)
    if file_name[-2:] != '.p' and file_name[-4:] != '.pkl':
        raise RuntimeError('Please use .p or .pkl for file extension')
    create_directory(file_name)
    with open(file_name, 'wb') as output_file:
        pickle.dump(image_obj, output_file, -1)

def load_image_object(object_file, image_file=None):
    """Load a pickled ImageAnalysis object."""
    if object_file[-2:] != '.p' and object_file[-4:] != '.pkl':
        raise RuntimeError('Please use .p or .pkl for file extension')
    object_file = true_path(object_file)
    with open(object_file, 'rb') as input_file:
        image_obj = pickle.load(input_file)
    if image_file is not None:
        image_file = true_path(image_file)
        image_obj.set_image_file(image_file)
    return image_obj

def create_directory(file_name):
    """Recursively creates directories if they don't exist."""
    file_name = true_path(file_name)
    file_list = file_name.split(os.sep)
    if '.' in file_list[-1]:
        file_list = file_list[:-1]
    directory = os.sep.join(file_list)
    if not os.path.exists(directory):
        print('Making directory', directory)
        os.makedirs(directory)

def save_data(image_obj, file_name):
    """Save object data as a dictionary in a text file."""
    file_name = true_path(file_name)
    data = to_dict(image_obj)
    if file_name[-3:] == 'txt':
        create_directory(file_name)
        with open(file_name, 'w') as save_file:
            save_file.write(str(data))
    else:
        raise RuntimeError('Please use .txt for file extension')

def load_data(file_name):
    """Return data from text file."""
    if file_name.endswith('.txt'):
        with open(file_name, 'r') as load_file:
            data = literal_eval(load_file.read())
    else:
        raise RuntimeError('Incorrect file type to load into object')
    return data

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

def true_path(path):
    """Use os to return the full path to a file."""
    return os.path.realpath(path)

def change_path(new_path, old_path, old_file):
    """Return the actual path based on the old relative path."""
    rel_path = os.path.relpath(old_file, old_path)
    return true_path(join_paths(new_path, rel_path))

def join_paths(path, *paths):
    """Join paths using os."""
    return os.sep.join((path,) + paths)

def get_directory(file_name=None):
    """Return the relevant directory."""
    if file_name is None:
        return os.getcwd()
    return os.path.split(true_path(file_name))[0] + os.sep
