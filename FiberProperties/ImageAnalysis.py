"""ImageAnalysis.py was written by Ryan Petersburg for use with fiber
characterization on the EXtreme PREcision Spectrograph
"""
import numpy as np
import cPickle as pickle
from ast import literal_eval
import os
from NumpyArrayHandler import *
from Calibration import Calibration

def getImageData(image_obj, method=None):
    """Returns relevant information from an ImageAnalysis object
    
    Args:
        image_obj [ImageAnalysis]: image object to be analyzed

    Returns:
        image_array [ndarray]: the 2D image
        y0 [float]: the fiber center y in pixels
        x0 [float]: the fiber center x in pixels
        radius [float]: the fiber radius in pixels
    """
    if method is None:
        if image_obj.getCamera() == 'ff':
            method = 'gaussian'
        elif image_obj.getCamera() == 'nf':
            method = 'edge'
        elif image_obj.getCamera() == 'in':
            method = 'edge'
        else:
            raise RuntimeError('Method or camera type must be declared to get image data')
    y0, x0 = image_obj.getFiberCenter(method=method, tol=1,
                                      test_range=10, show_image=False)
    radius = image_obj.getFiberRadius(method=method)        
    image_array = image_obj.getImage()
    return image_array, y0, x0, radius

def saveImageObject(image_obj, file_name):
    with open(file_name, 'wb') as output_file:
        pickle.dump(image_obj, output_file, -1)

def loadImageObject(file_name):
    with open(file_name, 'rb') as input_file:
        image_obj = pickle.load(input_file)
    return image_obj

class ImageAnalysis(object):
    """Fiber face image analysis class

    Class that conducts image analysis on a fiber face image after it has been
    corrected by the given dark, flat field, and ambient images. Also contains
    information about the CCD and camera that took the image. Public methods in
    this class allow calculation of the fiber face's centroid, center, and
    diameter using multiple different algorithms

    Args
    ----
    image_input : {None, 1D iterable, 2D iterable, string}
        The input used to set the image array. See
        NumpyArrayHandler.convertImageToArray() for details
    calibration : Calibration (default=None), optional
        Calibration object that contains the relevant dark, ambient, and flat
        fielded images in order to calibrate the analyzed image. If None,
        creates a Calibration object with all empty calibration images
    image_data : string, optional
        File location of previously calculated image data. Must be either a
        python pickle or text file containing a dictionary with image
        information formatted like the attributes in ImageAnalysis
    pixel_size : number, optional
        The side length of each CCD pixel in microns. This value should be
        contained in the image header, but if not, it needs to be defined
        upon initialization
    camera : {None, 'in', 'nf', 'ff'}, optional
        A string denoting the FCS camera type. Decides the magnification
        parameter (since it is known for the input and near field cameras),
        but otherwise is only used if printing or saving this information
        is necessary
    magnification : number, optional
        The magnification of the camera. Can also be set by choosing 'in' or
        'nf' for camera
    threshold : int, optional (default=256)
        The threshold value used with the centering method. This value may need
        to be tweaked for noisier images. Make sure to confirm decent results
        before setting this value too low
    kernel_size : odd int, optional (default=9)
        The kernel side length used when filtering the image. This value may
        need to be tweaked, especially with few co-added image, due to random
        noise. The filtered image is used for the centering algorithms, so for
        a "true test" using kernel_size=1, but be careful, because this may
        lead to needing a fairly high threshold for the noise.

    Attributes
    ----------
    image : 2D numpy.ndarray
    _calibration : Calibration
    _uncorrected_image : 2D numpy.ndarray
    _filtered_image : 2D numpy.ndarray
    
    _image_info : dict
    _analysis_info : dict

    _edges : dict
    _center : dict
    _centroid : dict
    _diameter : dict
    _array_sum : dict
    _fit : dict

    _phi : float

    """
    def __init__(self, image_input, calibration=None, image_data=None,
                 pixel_size=None, camera=None, magnification=None,
                 threshold=256, kernel_size=9):
        # Private attribute initialization 
        if image_data is None:
            self._image_info = dict(pixel_size=pixel_size,
                                    camera=camera,
                                    magnification=magnification,
                                    height=None,
                                    width=None,
                                    subframe_x=None,
                                    subframe_y=None,
                                    exp_time=None,
                                    bit_depth=None,
                                    date_time=None,
                                    temp=None,
                                    num_images=None,
                                    folder=None,
                                    test=None)
            self._analysis_info = dict(kernel_size=kernel_size,
                                       threshold=threshold)            
            self._edges = dict(left=None,
                               right=None,
                               top=None,
                               bottom=None)
            self._center = dict(edge=dict(x=None, y=None),
                                radius=dict(x=None, y=None),
                                circle=dict(x=None, y=None),
                                gaussian=dict(x=None, y=None))
            self._centroid = dict(edge=dict(x=None, y=None),
                                  radius=dict(x=None, y=None),
                                  circle=dict(x=None, y=None),
                                  gaussian=dict(x=None, y=None),
                                  full=dict(x=None, y=None))
            self._diameter = dict(edge=None,
                                  radius=None,
                                  circle=None,
                                  gaussian=None)
            self._array_sum = dict(radius=None,
                                   circle=None)
        else:
            self.loadData(image_data)

        self._fit = dict(gaussian=None,
                         polynomial=None)

        # Golden Ratio for optimization tests
        self._phi = (5 ** 0.5 - 1) / 2

        self._uncorrected_image = None
        self.image = None
        self._filtered_image = None        
        self._calibration = calibration

        if self._calibration is None:
            self._calibration = Calibration(None, None, None)
        self.setImageArray(image_input)

        self._filtered_image = self.getFilteredImage()

#=============================================================================#
#==== Private Variable Setters ===============================================#
#=============================================================================#

    def setImageArray(self, image_input):
        """Sets image to be analyzed

        Args
        ----
        image_input : {None, 1D iterable, 2D iterable, string}
            The input used to set the image array. See
            NumpyArrayHandler.convertImageToArray() for details

        Sets
        ----
        image

        """
        self._uncorrected_image, output_dict = convertImageToArray(image_input, True)
        if self._uncorrected_image is None:
            return

        self.setImageInfo(output_dict)

        self.image = self._calibration.executeErrorCorrections(self._uncorrected_image,
                                                               self._image_info['subframe_x'],
                                                               self._image_info['subframe_y'],
                                                               self._image_info['exp_time'])

    def setImageInfo(self, output_dict):
        """Sets image info using the output of convertImageToArray()
        
        Args
        ----
        output_dict : dict
            The dictionary of items to add to _image_info

        Sets
        ----
        _image_info[key]
            for all keys in output_dict.keys()
        """
        for key in output_dict:
            if self._image_info[key] is None:
                self._image_info[key] = output_dict[key]

        if self._image_info['magnification'] is None:
            if self._image_info['camera'] == 'nf' or self._image_info['camera'] == 'in':
                self._image_info['magnification'] = 10.0
            else:
                self._image_info['magnification'] = 1.0

#=============================================================================#
#==== Saving and Loading Data to File ========================================#
#=============================================================================#

    def loadData(self, file_name):
        """Loads data from a text file containing a python dictionary

        Args
        ----
        file_name : string
            The file where the data is located

        Raises
        ------
        RuntimeError
            if the file name does not end in '.p' or '.txt'

        Sets
        ----
        _image_info : dict
        _analysis_info : dict
        _edges : dict
        _center : dict
        _diameter : dict
        _centroid : dict
        _array_sum : dict

        """
        if file_name.endswith('.p'):
            data = pickle.load(open(file_name, 'rb'))
        elif file_name.endswith('.txt'):
            with open(file_name, 'r') as file:
                data = literal_eval(file.read())
        else:
            raise RuntimeError('Incorrect data type to load into object')

        self._image_info = data['image_info']
        self._analysis_info = data['analysis_info']
        self._edges = data['edges']
        self._center = data['center']
        self._diameter = data['diameter']
        self._centroid = data['centroid']
        self._array_sum = data['array_sum']

    def saveData(self, file_name=None, folder=None):
        """Pickle the data and also save the data as a text file dictionary

        Args
        ----
        folder : {None, string}, optional
            The containing folder to save the data. If None, uses the folder
            containing the first image initialized in the object
        file_name : {None, string}, optional
            The file name which is used to store the images. DO NOT include a
            file extension, as all files will be appended by '_data.txt' and
            '_data.p'

        Saves
        -----
        _image_info : dict
        _analysis_info : dict
        _edges : dict
        _center : dict
        _diameter : dict
        _centroid : dict
        _array_sum : dict

        """
        if file_name is None:
            file_name = self.getCamera()
        if folder is None:
            folder = self._image_info['folder']

        data = dict(image_info=self._image_info,
                    analysis_info=self._analysis_info,
                    edges=self._edges,
                    center=self._center,
                    diameter=self._diameter,
                    centroid=self._centroid,
                    array_sum=self._array_sum)

        file_base = folder + file_name
        
        pickle.dump(data, open(file_base + '_data.p', 'wb'))

        with open(file_base + '_data.txt', 'w') as file:
            file.write(str(data))

    def saveImages(self, file_name=None, folder=None):
        """Save image, uncorrected image, and filtered image as FITS images
        
        Args
        ----
        folder : {None, string}, optional
            The containing folder to save the images. If None, uses the folder
            containing the first image initialized in the object
        file_name : {None, string}, optional
            The file name which is used to store the images. DO NOT include a
            file extension, as all files will be saved as FITS and appended by
            the image type ('_uncorrected', '_corrected', and '_filtered')

        Saves
        -----
        image : 2D numpy.ndarray
            as FITS
        _uncorrected_image : 2D numpy.ndarray
            as FITS
        _filtered_image : 2D numpy.ndarray
            as FITS

        """
        if file_name is None:
            file_name = self.getCamera()
        if folder is None:
            folder = self._image_info['folder']

        file_base = folder + file_name

        #saveArray(self._uncorrected_image, file_base + '_uncorrected.fit')
        saveArray(self.image, file_base + '_corrected.fit')
        #saveArray(self._filtered_image, file_base + '_filtered.fit')

    def createDirectory(self, file_name):
        file_list = file_name.split('/')

        for i in xrange(len(file_list) - 2):
            if file_list[i + 1] not in os.listdir('/'.join(file_list[:i+1])):
                print file_list[i+1], file_list[:i+1]
                os.mkdir('/'.join(file_list[:i+2]))

#=============================================================================#
#==== Private Variable Getters ===============================================#
#=============================================================================#

    def convertPixelsToMicrons(self, value):
        """Converts a number or iterable from pixels to microns"""
        return convertPixelsToMicrons(value,
                                      self.getPixelSize(),
                                      self.getMagnification())

    def getFiberData(self, method=None, units='microns', **kwargs):
        """Return the fiber center and diameter

        Args
        ----
        method : {None, 'radius', 'gaussian', 'circle', 'edge'}, optional
            The method which is used to calculate the fiber center. If None,
            return the best calculated fiber center in the order 'radius' >
            'gaussian' > 'circle' > 'edge'
        units : {'pixels', 'microns'}, optional
            The units of the returned values
        **kwargs : 
            The keyworded arguments to pass to the centering method

        Sets
        ----
        _diameter[method] : float
            The diameter of the fiber face in the context of the given method
        _center[method] : {'x': float, 'y': float}
            The center of the fiber face in the context of the given method

        Returns
        -------
        _center[method]['y'] : float
            in the given units
        _center[method]['x'] : float
            in the given units
        _diameter[method] : float
            in the given units

        """
        center_y, center_x = self.getFiberCenter(method, units=units, **kwargs)
        kwargs['show_image'] = False
        diameter = self.getFiberDiameter(method, units=units, **kwargs)
        return center_y, center_x, diameter

    def getFiberRadius(self, method=None, units='pixels', **kwargs):
        """Return the fiber radius

        Finds the radius of the fiber using the given method or, if no method
        is given, the most precise method already completed

        Args
        ----
        method : {None, 'radius', 'gaussian', 'circle', 'edge'}, optional
            The method which is used to calculate the fiber center. If None,
            return the best calculated fiber center in the order 'radius' >
            'gaussian' > 'circle' > 'edge'
        units : {'pixels', 'microns'}, optional
            The units of the returned values
        **kwargs : 
            The keyworded arguments to pass to the centering method

        Sets
        ----
        _diameter[method] : float
            The diameter of the fiber face in the context of the given method
        _center[method] : {'x': float, 'y': float}
            The center of the fiber face in the context of the given method

        Returns
        -------
        _diameter[method] / 2.0 : float
            in the given units

        """
        return self.getFiberDiameter(method, units=units, **kwargs) / 2.0

    def getFiberDiameter(self, method=None, units='pixels', **kwargs):
        """Return the fiber diameter using the given method in the given units

        Find the diameter of the fiber using the given method or, if no method
        is given, the most precise method already completed

        Args
        ----
        method : {None, 'radius', 'gaussian', 'circle', 'edge'}, optional
            The method which is used to calculate the fiber center. If None,
            return the best calculated fiber center in the order 'radius' >
            'gaussian' > 'circle' > 'edge'
        units : {'pixels', 'microns'}, optional
            The units of the returned values
        **kwargs : 
            The keyworded arguments to pass to the centering method

        Sets
        ----
        _diameter[method] : float
            The diameter of the fiber face in the context of the given method
        _center[method] : {'x': float, 'y': float}
            The center of the fiber face in the context of the given method

        Returns
        -------
        _diameter[method] : float
            in the given units

        """
        if method is None:
            if self._diameter['radius'] is not None:
                method = 'radius'
            elif self._diameter['gaussian'] is not None:
                method = 'gaussian'
            else:
                method = 'edge'

        if self._diameter[method] is None:
            self.setFiberDiameter(method, **kwargs)

        diameter = self._diameter[method]

        if units == 'pixels':
            return diameter
        elif units == 'microns':
            return self.convertPixelsToMicrons(diameter)
        else:
            raise RuntimeError('Incorrect string for units')

    def getFiberCenter(self, method=None, units='pixels', **kwargs):
        """Return the fiber center using the given method in the given units

        Find the center position of the fiber using the given method or, if no
        method is given, the most precise method already completed

        Args
        ----
        method : {None, 'radius', 'gaussian', 'circle', 'edge'}, optional
            The method which is used to calculate the fiber center. If None,
            return the best calculated fiber center in the order 'radius' >
            'gaussian' > 'circle' > 'edge'
        units : {'pixels', 'microns'}, optional
            The units of the returned values
        **kwargs : 
            The keyworded arguments to pass to the centering method

        Sets
        ----
        _diameter[method] : float
            The diameter of the fiber face in the context of the given method
        _center[method] : {'x': float, 'y': float}
            The center of the fiber face in the context of the given method

        Returns
        -------
        _center[method]['y'] : float
            in the given units
        _center[method]['x'] : float
            in the given units

        """
        if method is None:
            if self._center['radius']['x'] is not None:
                method = 'radius'
            elif self._center['gaussian']['x'] is not None:
                method = 'gaussian'
            elif self._center['circle']['x'] is not None:
                method = 'circle'
            else:
                method = 'edge'

        if self._center[method]['x'] is None or method == 'circle':
            self.setFiberCenter(method, **kwargs)

        center = self._center[method]['y'], self._center[method]['x']

        if units == 'pixels':
            return center
        elif units == 'microns':
            return self.convertPixelsToMicrons(center)
        else:
            raise RuntimeError('Incorrect string for units')

    def getFiberCentroid(self, method=None, units='pixels', **kwargs):
        """Getter for the fiber centroid

        Args
        ----
        method : {None, 'full', 'edge', 'radius', 'gaussian', 'circle'}, optional
            See setFiberCentroid() for method details. If no method is given,
            chooses the most precise method already calculated in the order
            'radius' > 'gaussian' > 'circle' > 'edge' > 'full'
        units : {'pixels', 'microns'}, optional
            The units of the returned values
        **kwargs : 
            The keyworded arguments to pass to the centering and centroiding
            methods

        Sets
        ----
        _centroid[method] : {'x': float, 'y': float}
            The centroid of the fiber face in the context of the given method

        Returns
        -------
        _centroid[method]['y'] : float
            in the given units
        _centroid[method]['x'] : float
            in the given units
        """
        if method is None:
            if self._centroid['radius']['x'] is not None:
                method = 'radius'
            elif self._centroid['gaussian']['x'] is not None:
                method = 'gaussian'
            elif self._centroid['circle']['x'] is not None:
                method = 'circle'
            elif self._centroid['edge']['x'] is not None:
                method = 'edge'
            else:
                method = 'full'

        if self._centroid[method]['x'] is None:
            self.setFiberCentroid(method, **kwargs)
        centroid = (self._centroid[method]['y'], self._centroid[method]['x'])

        if units == 'pixels':
            return centroid
        elif units == 'microns':
            return self.convertPixelsToMicrons(centroid)
        else:
            raise RuntimeError('Incorrect string for units')

    def getGaussianFit(self):
        """Return the best gaussian fit for the image

        Returns
        -------
        _fit['gaussian'] : 2D numpy.ndarray

        """
        if self._fit['gaussian'] is None:
            self.setFiberCenter(method='gaussian')
        return self._fit['gaussian']

    def getMeshGrid(self):
        """Return a meshgrid of the same size as the stored image"""
        return meshGridFromArray(self.image)

    def getPolynomialFit(self, deg=6, x0=None, y0=None):
        """Return the best polynomial fit for the image
        
        Args
        ----
        deg : int (default=6)
            The degree of polynomial to fit
        x0 : number
            The center column to use for the radial polynomial. Uses best
            calculated center if None.
        y0 : number
            The center row to use for the radial polynomial. Uses best
            calculated center if None.

        Returns
        -------
        polynomial_fit : 2D numpy.ndarray

        """
        if self._fit['polynomial'] is None:
            if y0 is None or x0 is None:
                y0, x0 = self.getFiberCenter()
            self._fit['polynomial'] = polynomialFit(self.image, deg, x0, y0)
        return self._fit['polynomial']

    def getTophatFit(self):
        """Return the circle array that best covers the fiber face
        
        Returns
        -------
        circle_array : numpy.ndarray (2D)
            Circle array centered at best calculated center and with best
            calculated diameter

        """
        y0, x0 = self.getFiberCenter(show_image=False)
        radius = self.getFiberRadius(show_image=False)
        return self.image.max() * circleArray(self.getMeshGrid(), x0, y0, radius, res=1)

    def getFilteredImage(self, kernel_size=None):
        """Return a median filtered image
        
        Args
        ----
        kernel_size : {None, int (odd)}, optional
            The side length of the kernel used to median filter the image. Uses
            self._analysis_info['kernel_size'] if None.
        
        Returns
        -------
        filtered_image : 2D numpy.ndarray
            The stored image median filtered with the given kernel_size
        None : NoneType
            If self.image is None

        """
        if self.image is None:
            return None
        if kernel_size is None and self._filtered_image is not None:
            return self._filtered_image
        if kernel_size is None:
            kernel_size = self._analysis_info['kernel_size']
        return filteredImage(self.image, kernel_size)

    def getDarkImage(self):
        """Return the dark image from the calibration object"""
        return self._calibration.dark_image

    def getAmbientImage(self):
        """Return the ambient image from the calibration object"""
        return self._calibration.ambient_image

    def getFlatImage(self):
        """Return the flatfield image from the calibration object"""
        return self._calibration.flat_image

    def getArraySum(self):
        """Return the sum of the stored image array"""
        return sumArray(image_array)

    def getImageInfo(self, info_type=None):
        """Return the image info dictionary or contained quantity

        Args
        ----
        info_type : {None, string in _analysis_info.keys()}, optional
            If None, return the entire _image_info dictionary. Otherwise return
            the information under _analysis_info[info_type]

        Returns
        -------
        _image_info : dict
            When info_type not given or None
        _image_info[info_type] : {number, string, dict}

        Raises
        ------
        RuntimeError
            if an incorrect info_type is given

        """
        if info_type is None:
            return self._image_info
        elif info_type in self._image_info:
            if self._image_info[info_type] is None:
                raise RuntimeError(info_type + ' needs to be set externally')
            return self._image_info[info_type]
        raise RuntimeError('Incorrect string for image info property')

    def getAnalysisInfo(self, info_type=None):
        """Return the analysis info dictionary or contained quantity
        
        Args
        ----
        info_type : {None, string in _analysis_info.keys()}, optional
            If None, return the entire _analysis_info dictionary. Otherwise
            return the information under _analysis_info[info_type]

        Returns
        -------
        _analysis_info : dict
            When info_type is not given or None
        _analysis_info[info_type] : number

        Raises
        ------
        RuntimeError
            if an incorrect info_type is given

        """
        if info_type is None:
            return self._analysis_info
        elif info_type in self._analysis_info:
            if self._analysis_info[info_type] is None:
                raise RuntimeError(info_type + ' needs to be set externally')
            return self._analysis_info[info_type]
        raise RuntimeError('Incorrect string for image info property')

    def getHeight(self):
        """Return the image height in pixels"""
        return self._image_info['height']

    def getWidth(self):
        """Return the image width in pixels"""
        return self._image_info['width']

    def getMagnification(self):
        """Return the magnification"""
        if self._image_info['magnification'] is None:
            raise RuntimeError('Magnification needs to be set externally')
        return self._image_info['magnification']

    def getCamera(self):
        """Return the string denoting the camera type"""
        if self._image_info['camera'] is None:
            raise RuntimeError('Camera needs to be set externally')
        return self._image_info['camera']

    def getPixelSize(self):
        """Return the pixel size in microns"""
        if self._image_info['pixel_size'] is None:
            raise RuntimeError('Pixel Size needs to be set externally')
        return self._image_info['pixel_size']

    def getNumImages(self):
        """Return the number of images that created the object"""
        if self._image_info['num_images'] is None:
            raise RuntimeError('Number of images improperly set')
        return self._image_info['num_images']

    def getImage(self):
        """Return the image"""
        if self.image is None:
            raise RuntimeError('Images were not set for this object')
        return self.image

#=============================================================================#
#==== Image Centroiding ======================================================#
#=============================================================================#

    def setFiberCentroid(self, method='full', radius_factor=None,
                         show_image=False, **kwargs):
        """Find the centroid of the fiber face image

        Args
        ----
        method : {'full', 'edge', 'radius', 'gaussian', 'circle'}, optional
            If 'full', takes the centroid of the entire image. Otherwise, uses
            the specified method to isolate only the fiber face in the image
        radius_factor : number
            The factor by which the radius is multiplied when isolating the
            fiber face in the image

        Sets
        ----
        _centroid[method] : {'x': float, 'y': float}
            The centroid of the image in the context of the given method

        """
        if method != 'full':
            if radius_factor is None:
                radius_factor = 1.0

            y0, x0 = self.getFiberCenter(method=method, show_image=False, **kwargs)
            radius = self.getFiberRadius(method=method, show_image=False, **kwargs)
            image_array_iso = isolateCircle(self.image, x0, y0,
                                            radius*radius_factor, res=1)

        else:
            image_array_iso = self.image

        x_array, y_array = self.getMeshGrid()
        self._centroid[method]['x'] = (image_array_iso * x_array).sum() / image_array_iso.sum()
        self._centroid[method]['y'] = (image_array_iso * y_array).sum() / image_array_iso.sum()

        if show_image:
            plotDot(self.image, self._centroid[method]['y'], self._centroid[method]['x'])

#=============================================================================#
#==== Image Centering ========================================================#
#=============================================================================#

    def setFiberData(self, method, **kwargs):
        """Set the fiber center, diameter, and centroid using the same method

        Args
        ----
        method : {'edge', 'radius', 'gaussian', 'circle'}
            Uses the respective method to find the fiber center
        **kwargs : 
            The keyworded arguments to pass to the centering method

        Sets
        ----
        _centroid[method] : {'x': float, 'y': float}
            The centroid of the image in the context of the given method
        _center[method] : {'x': float, 'y': float}
            The center of the fiber face in the context of the given method
        _diameter[method] : float
            The diameter of the fiber face in the context of the given method

        """
        self.setFiberCenter(method, **kwargs)
        self.setFiberCentroid(method, **kwargs)

    def setFiberDiameter(self, method, **kwargs):
        """Set the fiber diameter using given method

        Args
        ----
        method : {'edge', 'radius', 'gaussian', 'circle'}
            Uses the respective method to find the fiber center
        **kwargs :
            The keyworded arguments to pass to the centering method

        Sets
        ----
        _diameter[method] : float
            The diameter of the fiber face in the context of the given method
        _center[method] : {'x': float, 'y': float}
            The center of the fiber face in the context of the given method

        Raises
        ------
        RuntimeError
            cannot accept the 'circle' method when setting the diameter since
            it requires a known radius to run
        """
        if method == 'circle':
            raise RuntimeError('Fiber diameter cannot be set by circle method')
        self.setFiberCenter(method, **kwargs)

    def setFiberCenter(self, method, show_image=False, **kwargs):
        """Find fiber center using given method

        Args
        ----      
        method : {'edge', 'radius', 'gaussian', 'circle'}
            Uses the respective method to find the fiber center
        show_image : boolean, optional (default=False)
            Whether or not to show relevant fitting image
        **kwargs :
            The keyworded arguments to pass to the centering method

        Raises
        ------
        RuntimeError
            needs a valid method string to run the proper algorithm
        """
        # Reset the fits due to new fiber parameters
        self._fit['gaussian'] = None
        self._fit['polynomial'] = None

        if method == 'radius':
            self.setFiberCenterRadiusMethod(**kwargs)
        elif method == 'edge':
            self.setFiberCenterEdgeMethod()
        elif method == 'circle':
            self.setFiberCenterCircleMethod(**kwargs)
        elif method == 'gaussian':
            self.setFiberCenterGaussianMethod()
        else:
            raise RuntimeError('Incorrect string for fiber centering method')

        if show_image:
            if method == 'gaussian':
                showImageArray(self._fit['gaussian'])
                plotOverlaidCrossSections(self.image, self.getGaussianFit(),
                                          self._center['gaussian']['y'], 
                                          self._center['gaussian']['x'])
            else:
                self.showOverlaidTophat(self._center[method]['x'],
                                        self._center[method]['y'],
                                        self._diameter[method] / 2.0,
                                        tol=1)

    def setFiberCenterGaussianMethod(self):
        """Set fiber center using a Gaussian Fit

        Uses Scipy.optimize.curve_fit method to fit fiber image to
        gaussianArray(). The radius found extends to 2-sigma of the gaussian
        therefore encompassing ~95% of the imaged light. Use previous methods
        of center-finding to approximate the location of the center

        Sets
        ----
        _diameter['gaussian'] : float
            Diameter of the fiber in the gaussian method context
        _center['gaussian'] : {'x': float, 'y': float}
            Center of the fiber in the gaussian method context
        _fit['gaussian'] : 2D numpy.ndarray
            Best gaussian fit for the fiber image

        """
        y0, x0 = self.getFiberCenter(method='edge')
        radius = self.getFiberRadius(method='edge')

        if self._image_info['camera'] == 'in':
            approx_circle_array = (np.median(intensityArray(self._filtered_image, x0, y0, radius))
                                   * circleArray(self.getMeshGrid(), x0, y0, radius))

            filtered_image, new_x0, new_y0 = cropImage(self._filtered_image - approx_circle_array,
                                                       x0, y0, radius)

            initial_guess = (new_x0, new_y0, 100 / self.getPixelSize(),
                             filtered_image.max(), filtered_image.min())
            _, opt_parameters = gaussianFit(filtered_image,
                                            initial_guess=initial_guess,
                                            full_output=True)

            x0 = opt_parameters[0] + int(x0-radius)
            y0 = opt_parameters[1] + int(y0-radius)
            radius = abs(opt_parameters[2])
            amp = opt_parameters[3]
            offset = opt_parameters[4]

            self._fit['gaussian'] = (gaussianArray(self.getMeshGrid(), x0, y0,
                                                  radius, amp,
                                                  offset).reshape(self.getHeight(),
                                                                  self.getWidth())
                                     + approx_circle_array)

        else:
            initial_guess = (x0, y0, radius, self._filtered_image.max(),
                             self._filtered_image.min())

            self._fit['gaussian'], opt_parameters = gaussianFit(self._filtered_image,
                                                                initial_guess=initial_guess,
                                                                full_output=True)
            x0 = opt_parameters[0]
            y0 = opt_parameters[1]
            radius = abs(opt_parameters[2])

        self._center['gaussian']['x'] = x0
        self._center['gaussian']['y'] = y0
        self._diameter['gaussian'] = radius * 2.0

    def setFiberCenterRadiusMethod(self, tol=1, test_range=None):
        """Set fiber center using dark circle with varying radius

        Uses a golden mean optimization method to find the optimal radius of the
        dark circle that covers the fiber image used in
        getFiberCenterCircleMethod(). The optimization is for a parameter
        array_sum which is weighted by the area of the circle, meaning that a
        smaller circle is preferred over one that simply covers the entire image

        Args
        ----
        tol : number (default=1)
            Minimum possible range of radius values before ending iteration
        test_range: int (in pixels)
            Range of tested radii, i.e. max(radius) - min(radius). If None,
            uses full possible range

        Sets
        ----
        _diameter['radius'] : float
            Diameter of the fiber in the radius method context
        _center['radius'] : {'x': float, 'y': float}
            Center of the fiber in the radius method context
        _diameter['circle'] : float
            Also uses the circle method, therefore changes this value
        _center['circle'] : float
            Also uses the circle method, therefore chnages this value 

        """
        # Initialize range of tested radii
        r = np.zeros(4).astype(float)

        if test_range is not None:
            approx_radius = self.getFiberRadius()
            test_range = test_range / 2.0

            r[0] = approx_radius - test_range
            if r[0] < 0.0:
                r[0] = 0.0
            r[3] = approx_radius + test_range
        else:
            r[0] = 0
            r[3] = min(self._image_info['height'], self._image_info['width']) / 2.0

        r[1] = r[0] + (1 - self._phi) * (r[3] - r[0])
        r[2] = r[0] + self._phi * (r[3] - r[0])

        array_sum = np.zeros(2).astype(float)
        for i in xrange(2):
            self.setFiberCenter(method='circle', radius=r[i+1],
                                tol=tol, test_range=test_range)
            array_sum[i] = (self._array_sum['circle']
                            + self._analysis_info['threshold']
                            * np.pi * r[i+1]**2)

        min_index = np.argmin(array_sum) # Integer 0 or 1 for min of r[1], r[2]

        while abs(r[3]-r[0]) > tol:
            if min_index == 0:
                r[3] = r[2]
                r[2] = r[1]
                r[1] = r[0] + (1 - self._phi) * (r[3] - r[0])
            else:
                r[0] = r[1]
                r[1] = r[2]
                r[2] = r[0] + self._phi * (r[3] - r[0])

            array_sum[1 - min_index] = array_sum[min_index]

            self.setFiberCenter(method='circle', radius=r[min_index+1],
                                tol=tol, test_range=test_range)
            array_sum[min_index] = (self._array_sum['circle']
                                    + self._analysis_info['threshold']
                                    * np.pi * r[min_index+1]**2)

            min_index = np.argmin(array_sum) # Integer 0 or 1 for min of r[1], r[2]

        self._diameter['radius'] = r[min_index+1] * 2
        self._center['radius']['y'] = self._center['circle']['y']
        self._center['radius']['x'] = self._center['circle']['x']
        self._array_sum['radius'] = np.amin(array_sum)

    def setFiberCenterCircleMethod(self, radius, tol=1, test_range=None):
        """Finds fiber center using a dark circle of set radius

        Uses golden mean method to find the optimal center for a circle
        covering the fiber image. The optimization is for a parameter array_sum
        that simply sums over the entire fiber image array

        Args
        ----
        radius : float
            Radius to use when creating circle
        tol : number (default=1)
            Minimum possible range of center values before ending iteration
        test_range: int (in pixels)
            Range of tested centers, i.e. max(x0) - min(x0). If None,
            uses full possible range

        Sets
        ----
        _diameter['circle'] : float
            Diameter of the fiber in the circle method context
        _center['circle'] : {'x': float, 'y': float}
            Center of the fiber in the circle method context
        _diameter['edge'] : float
            If test_range is not None, approximates the circle's center using
            the edge method
        _center['edge'] : float
            If test_range is not None, approximates the circle's center using
            the edge method

        """
        #print 'Testing diameter '
        res = int(1.0/tol)

        # Create four "corners" to test center of the removed circle
        x = np.zeros(4).astype(float)
        y = np.zeros(4).astype(float)

        if test_range is not None:
            approx_center = self.getFiberCenter(method='edge', show_image=False)
            test_range = test_range / 2.0

            x[0] = approx_center[1] - test_range
            if x[0] < radius:
                x[0] = radius
            x[3] = approx_center[1] + test_range
            if x[3] > self._image_info['width'] - radius:
                x[3] = self._image_info['width'] - radius

            y[0] = approx_center[0] - test_range
            if y[0] < radius:
                y[0] = radius
            y[3] = approx_center[0] + test_range
            if y[3] > self._image_info['height'] - radius:
                y[3] = self._image_info['height'] - radius

        else:
            x[0] = radius
            x[3] = self._image_info['width'] - radius

            y[0] = radius
            y[3] = self._image_info['height'] - radius

        x[1] = x[0] + (1 - self._phi) * (x[3] - x[0])
        x[2] = x[0] + self._phi * (x[3] - x[0])

        y[1] = y[0] + (1 - self._phi) * (y[3] - y[0])
        y[2] = y[0] + self._phi * (y[3] - y[0])

        # Initialize array sums to each corner
        array_sum = np.zeros((2, 2)).astype(float)
        for i in xrange(2):
            for j in xrange(2):
                removed_circle_array = removeCircle(self._filtered_image,
                                                    x[i+1], y[j+1],
                                                    radius, res)
                array_sum[j, i] = sumArray(removed_circle_array)

        # Find the index of the corner with minimum array_sum
        min_index = np.unravel_index(np.argmin(array_sum), (2, 2)) # Tuple

        while abs(x[3] - x[0]) > tol and abs(y[3] - y[0]) > tol:
            # Move the other corners to smaller search area
            if min_index[0] == 0:
                y[3] = y[2]
                y[2] = y[1]
                y[1] = y[0] + (1 - self._phi) * (y[3] - y[0])
            else:
                y[0] = y[1]
                y[1] = y[2]
                y[2] = y[0] + self._phi * (y[3] - y[0])
            if min_index[1] == 0:
                x[3] = x[2]
                x[2] = x[1]
                x[1] = x[0] + (1 - self._phi) * (x[3] - x[0])
            else:
                x[0] = x[1]
                x[1] = x[2]
                x[2] = x[0] + self._phi * (x[3] - x[0])

            # Replace the opposite corner array sum (so it doesn't need to be recalculated)
            array_sum[1 - min_index[0], 1 - min_index[1]] = array_sum[min_index]
            min_index = (1 - min_index[0], 1 - min_index[1])

            # Recalculate new sums for all four corners
            for i in xrange(2):
                for j in xrange(2):
                    if i != min_index[1] or j != min_index[0]:
                        removed_circle_array = removeCircle(self._filtered_image,
                                                            x[i+1], y[j+1],
                                                            radius, res)
                        array_sum[j, i] = sumArray(removed_circle_array)

            min_index = np.unravel_index(np.argmin(array_sum), (2, 2))

        self._center['circle']['x'] = x[min_index[1]+1]
        self._center['circle']['y'] = y[min_index[0]+1]
        self._diameter['circle'] = radius * 2.0
        self._array_sum['circle'] = np.amin(array_sum)

    def setFiberCenterEdgeMethod(self):
        """The averages of the fiber edges gives the fiber center

        Returns:
            center_y, center_x
        """
        self.setFiberEdges()

        self._center['edge']['y'] = (self._edges['top'] + self._edges['bottom']) / 2.0
        self._center['edge']['x'] = (self._edges['left'] + self._edges['right']) / 2.0

    def setFiberEdges(self):
        """Set fiber edge pixel values

        Sets the left, right, top, and bottom edges of the fiber by finding where
        the maxima of each row and column cross the given threshold. Also sets
        the width of the fiber by the maximum of the horizontal and vertical
        lengths

        Sets:
            self._edges['left']
            self._edges['right']
            self._edges['top']
            self._edges['bottom']
            self._diameter['edge']
        """
        left = -1
        right = -1
        for index in xrange(self._image_info['width']):
            if left < 0:
                if self._filtered_image[:, index].max() > self._analysis_info['threshold']:
                    left = index
            else:
                if self._filtered_image[:, index].max() > self._analysis_info['threshold']:
                    right = index

        top = -1
        bottom = -1
        for index in xrange(self._image_info['height']):
            if top < 0:
                if self._filtered_image[index, :].max() > self._analysis_info['threshold']:
                    top = index
            else:
                if self._filtered_image[index, :].max() > self._analysis_info['threshold']:
                    bottom = index

        self._edges['left'] = left
        self._edges['right'] = right
        self._edges['top'] = top
        self._edges['bottom'] = bottom
        self._diameter['edge'] = ((right - left) + (bottom - top)) / 2.0

#=============================================================================#
#==== Overriding Methods =====================================================#
#=============================================================================#

    def plotCrossSections(self, image_array=None, row=None, column=None):
        if image_array is None:
            image_array = self.image
        if row is None:
            row = self._image_info['height'] / 2.0
        if column is None:
            column = self._image_info['width'] / 2.0
        plotCrossSections(image_array, row, column)

    def showImageArray(self, image_array=None):
        if image_array is None:
            image_array = self.image
        showImageArray(image_array)


    def showOverlaidTophat(self, x0, y0, radius, tol=1):
        res = int(1.0/tol)
        showImageArray(removeCircle(self.image, x0, y0, radius, res=res))
        plotOverlaidCrossSections(self.image,
                                  2 * self._analysis_info['threshold']
                                  *circleArray(self.getMeshGrid(),
                                               x0, y0, radius, res=res),
                                  y0, x0)


if __name__ == "__main__":
    folder = 'Stability Measurements/2016-08-15 Stability Test Unagitated/'

    calibration = Calibration(dark=[folder + 'Dark/in_' + str(i).zfill(3) + '.fit' for i in xrange(10)],
                              ambient=[folder + 'Ambient/in_' + str(i).zfill(3) + '.fit' for i in xrange(10)],
                              flat=[folder + 'Flat/in_' + str(i).zfill(3) + '.fit' for i in xrange(8)])

    images = [folder + 'Images/in_' + str(i).zfill(3) + '.fit' for i in xrange(100)]

    min_image = 0
    max_image = 100
    im_obj = ImageAnalysis(images[min_image:max_image+1], calibration, magnification=1)

    tol = 0.25
    test_range = 5
    factor = 1.0

    im_obj.showImageArray()
    for key in im_obj._image_info:
        print key + ': ' + str(im_obj._image_info[key])
    for key in im_obj._analysis_info:
        print key + ': ' + str(im_obj._analysis_info[key])
    print
    print 'Centroid'
    centroid_row, centroid_column = im_obj.getFiberCentroid(method='full', radius_factor=factor)
    print 'Centroid Row:', centroid_row, 'Centroid Column:', centroid_column
    print
    print 'Edge:'
    center_y, center_x = im_obj.getFiberCenter(method='edge', show_image=True)
    centroid_row, centroid_column = im_obj.getFiberCentroid(method='edge', radius_factor=factor)
    print 'Diameter:', im_obj.getFiberDiameter(method='edge', units='microns'), 'microns'
    print 'Center Row:', center_y, 'Center Column:', center_x
    print
    print 'Radius:'
    center_y, center_x = im_obj.getFiberCenter(method= 'radius', tol=tol, test_range=test_range, show_image=True)
    centroid_row, centroid_column = im_obj.getFiberCentroid(method='radius', radius_factor=factor)
    print 'Diameter:', im_obj.getFiberDiameter(method='radius', units='microns'), 'microns'
    print 'Center Row:', center_y, 'Center Column:', center_x
    print
    print 'Gaussian:'
    center_y, center_x = im_obj.getFiberCenter(method='gaussian', show_image=True)
    centroid_row, centroid_column = im_obj.getFiberCentroid(method='gaussian', radius_factor=factor)
    print 'Diameter:', im_obj.getFiberDiameter(method='gaussian', units='microns'), 'microns'
    print 'Center Row:', center_y, 'Center Column:', center_x

    im_obj.saveData(file_name='in', folder=folder)
    im_obj.saveImages(file_name='in', folder=folder)

    im_obj = ImageAnalysis(images, calibration, image_data=folder + 'in_data.p')

    for key in im_obj._image_info:
        print key + ': ' + str(im_obj._image_info[key])
    for key in im_obj._analysis_info:
        print key + ': ' + str(im_obj._analysis_info[key])
    print
    print 'Centroid'
    centroid_row, centroid_column = im_obj.getFiberCentroid(method='full', radius_factor=factor)
    print 'Centroid Row:', centroid_row, 'Centroid Column:', centroid_column
    print
    print 'Edge:'
    center_y, center_x = im_obj.getFiberCenter(method='edge')
    centroid_row, centroid_column = im_obj.getFiberCentroid(method='edge', radius_factor=factor)
    print 'Diameter:', im_obj.getFiberDiameter(method='edge', units='microns'), 'microns'
    print 'Center Row:', center_y, 'Center Column:', center_x
    print
    print 'Radius:'
    center_y, center_x = im_obj.getFiberCenter(method= 'radius', tol=tol, test_range=test_range)
    centroid_row, centroid_column = im_obj.getFiberCentroid(method='radius', radius_factor=factor)
    print 'Diameter:', im_obj.getFiberDiameter(method='radius', units='microns'), 'microns'
    print 'Center Row:', center_y, 'Center Column:', center_x
    print
    print 'Gaussian:'
    center_y, center_x = im_obj.getFiberCenter(method='gaussian')
    centroid_row, centroid_column = im_obj.getFiberCentroid(method='gaussian', radius_factor=factor)
    print 'Diameter:', im_obj.getFiberDiameter(method='gaussian', units='microns'), 'microns'
    print 'Center Row:', center_y, 'Center Column:', center_x