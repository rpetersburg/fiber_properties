"""ImageAnalysis.py was written by Ryan Petersburg for use with fiber
characterization on the EXtreme PREcision Spectrograph
"""
import numpy as np
import re
from PIL import Image
from astropy.io import fits
from datetime import datetime
from collections import Iterable
from Containers import ImageInfo

def convert_image_to_array(image_input, full_output=False):
    """Converts an image input to a numpy array or None

    Args
    ----
    image_input : {None, 1D iterable, 2D iterable, string}
        Inputting None simply returns None. Inputting a string of a file name
        returns the image contained within that file. Inputting an iterable
        containing strings returns all of the images in those files co-added
        together. Inputting a 2D iterable returns a 2D numpy.ndarray of the
        input iterable. Inputting a 1D iterable containing 2D iterables returns
        those 2D iterables co-added together in a single numpy.ndarray
    full_output : boolean, optional (default=False)
        Whether or not to include relevant information from the image header in
        the return. Automatically False if the image_input is an ndarray (and
        therefore without a header).

    Returns
    -------
    image_array : 2D numpy.ndarray or None
        2D numpy array if the image input checks out, None otherwise
    image_info : ImageInfo, optional
        ImageInfo object containing information from the image header
    """
    image_array = None

    if image_input is None:
        pass

    # Image input is a single file name
    elif isinstance(image_input, basestring):
        image_array = image_array_from_file(image_input, full_output)
        if full_output:
            image_info = image_array[1]
            image_array = image_array[0]
            image_info.num_images = 1

    # Image input is a sequence of file names
    elif isinstance(image_input, Iterable) and isinstance(image_input[0], basestring):
        list_len = float(len(image_input))
        image_array = image_array_from_file(image_input[0], full_output)
        if full_output:
            image_info = image_array[1]
            image_array = image_array[0]
            image_info.num_images = list_len
        image_array /= list_len
        for image_string in image_input[1:]:
            image_array += image_array_from_file(image_string) / list_len

    # Image input is a single array
    elif isinstance(image_input, Iterable) and len(np.array(image_input).shape) == 2:
        image_array = np.array(image_input)
        if full_output:
            image_info.num_images = 1

    # Image input is a sequence of arrays
    elif isinstance(image_input, Iterable) and isinstance(image_input[0], Iterable):
        list_len = float(len(image_input))
        image_input = np.array(image_input)
        image_array = image_input[0] / list_len
        for image in image_input[1:]:
            image_array += image / list_len            
        if full_output:
            image_info.num_images = list_len
        
    else:
        raise RuntimeError('Incorrect type for image input')

    if full_output:
        if image_array is not None:
            image_info.height, image_info.width = image_array.shape
            return image_array, image_info
        return image_array, ImageInfo()
    return image_array

def image_array_from_file(image_string, full_output=False):
    """Returns image from file as 2D np.ndarray
    
    Args
    ----
    image_string : string
        File location of the image (FITS or TIFF) to be converted
    full_output : boolean, optional
        whether or not to include relevant information from the image header in
        the return

    Returns
    -------
    image_array : 2D numpy.ndarray
        2D numpy array of the file's image
    image_info : ImageInfo
        Object containing the information from the image header

    """
    if image_string[-3:] == 'fit':
        image = fits.open(image_string)[0]
        image_array = image.data.astype('float64')
        if full_output: 
            header = dict(image.header)

    elif image_string[-3:] == 'tif':
        image = Image.open(image_string)
        image_array = np.array(image).astype('float64')
        if full_output:
            # Complicated way to get the header from a TIF image as a dictionary
            header = dict([i.split('=') for i in image.tag[270][0].split('\r\n')][:-1])
            header['BITPIX'] = int(image.tag[258][0])

    else:
        raise ValueError('Incorrect image file extension')

    if full_output:
        image_info = ImageInfo()   
        image_info.folder = '/'.join(image_string.split('/')[:-1]) + '/'

        if 'XORGSUBF' in header:
            image_info.subframe_x = int(header['XORGSUBF'])
            image_info.subframe_y = int(header['YORGSUBF'])

        image_info.bit_depth = int(header['BITPIX'])
        if 'XPIXSZ' in header:
            image_info.pixel_size = float(header['XPIXSZ'])
        if 'EXPTIME' in header:
            image_info.exp_time = float(header['EXPTIME'])
        if 'DATE-OBS' in header:
            try:
                image_info.date_time = datetime.strptime(header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')
            except ValueError:
                image_info.date_time = datetime.strptime(header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S')
        if 'CCD-TEMP' in header:
            image_info.temp = float(header['CCD-TEMP'])

        image_string_list = image_string.split('/')
        if 'TELESCOP' in header:
            image_info.camera = str(header['TELESCOP'])
        elif 'nf' in image_string_list or 'nf_' in image_string_list[-1]:
            image_info.camera = 'nf'
        elif 'ff' in image_string_list or 'ff_' in image_string_list[-1]:
            image_info.camera = 'ff'
        elif 'in' in image_string_list or 'in_' in image_string_list[-1]:
            image_info.camera = 'in'

        if 'OBJECT' in header:
            image_info.test = str(header['OBJECT'])

        return image_array, image_info

    return image_array
