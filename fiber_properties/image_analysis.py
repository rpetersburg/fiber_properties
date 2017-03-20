"""image_analysis.py was written by Ryan Petersburg for use with fiber
characterization on the EXtreme PREcision Spectrograph
"""
from ast import literal_eval
from collections import Iterable
import numpy as np

from fiber_properties.numpy_array_handler import (sum_array,
                                                  mesh_grid_from_array,
                                                  crop_image, remove_circle,
                                                  isolate_circle, filter_image,
                                                  gaussian_array, circle_array,
                                                  polynomial_fit, gaussian_fit,
                                                  rectangle_array, rectangle_fit)
from fiber_properties.plotting import (plot_cross_sections,
                                       plot_overlaid_cross_sections, plot_dot,
                                       show_image_array, show_plots,
                                       plot_image_array)
from fiber_properties.input_output import (save_array, save_image_object,
                                           save_data)
from fiber_properties.image_conversion import convert_image_to_array
from fiber_properties.calibration import Calibration
from fiber_properties.containers import (ImageInfo, AnalysisInfo, FiberInfo,
                                         Edges, FRDInfo)

#=============================================================================#
#===== Useful Functions ======================================================#
#=============================================================================#

def get_image_data(image_obj, method=None):
    """Returns relevant information from an ImageAnalysis object

    Args
    ----
    image_obj [ImageAnalysis]: image object to be analyzed

    Returns
    -------
    image_array : ndarray
        the 2D image
    y0 : float
        the fiber center y in pixels
    x0 : float
        the fiber center x in pixels
    radius : float
        the fiber radius in pixels
    """
    if method is None:
        if image_obj.get_camera() == 'ff':
            method = 'gaussian'
        elif image_obj.get_camera() == 'nf':
            method = 'edge'
        elif image_obj.get_camera() == 'in':
            method = 'edge'
        else:
            raise RuntimeError('Method or camera type must be declared to get image data')
    y0, x0 = image_obj.get_fiber_center(method=method, tol=1,
                                        test_range=10, show_image=False)
    radius = image_obj.get_fiber_radius(method=method)
    image_array = image_obj.get_image()
    return image_array, y0, x0, radius

def convert_pixels_to_microns(value, pixel_size, magnification):
    """Converts a value or iterable from pixels to microns"""
    if isinstance(value, Iterable):
        return tuple(np.array(value) * pixel_size / magnification)
    return value * pixel_size / magnification

def convert_pixels_to_units(value, pixel_size, magnification, units):
    """Converts a value or iterable from pixels to given units"""
    if units == 'pixels':
        return value
    elif units == 'microns':
        return convert_pixels_to_microns(value, pixel_size, magnification)
    else:
        raise RuntimeError('Incorrect string for units')

def convert_microns_to_units(value, pixel_size, magnification, units):
    """Converts a value or iterable from microns to given units"""
    if units == 'microns':
        return value
    elif units == 'pixels':
        if isinstance(value, Iterable):
            return tuple(np.array(value) * magnification / pixel_size)
        return value * magnification / pixel_size

def convert_fnum_to_radius(fnum, pixel_size, magnification, units='pixels'):
    """Converts a focal ratio to an image radius in given units."""
    fcs_focal_length = 4.0 # inches
    diameter = fcs_focal_length / fnum # inches
    radius = 25400 * diameter / 2.0 # microns
    return convert_microns_to_units(radius, pixel_size, magnification, units)

#=============================================================================#
#===== ImageAnalysis Class ===================================================#
#=============================================================================#

class ImageAnalysis(object):
    """Fiber face image analysis class

    Class that conducts image analysis on a fiber face image after it has been
    corrected by the given dark, flat field, and ambient images. Also contains
    information about the CCD and camera that took the image. Public methods in
    this class allow calculation of the fiber face's centroid, center, and
    diameter using multiple different algorithms

    Args
    ----
    image_input : {None, string, list(string), 2D numpy.ndarray,
                   list(2D numpy.ndarray)}
        The input used to set the image array. See
        ImageConversion.convert_image_to_array() for details
    dark : {None, string, list(string), 2D numpy.ndarray,
            list(2D numpy.ndarray)}, optional
        The input used to set the dark image. See Calibration module and
        ImageConcersion.convert_image_to_array() for details
    ambient : {None, string, list(string), 2D numpy.ndarray,
               list(2D numpy.ndarray)}, optional
        The input used to set the ambient image. See Calibration module and
        ImageConcersion.convert_image_to_array() for details
    flat : {None, string, list(string), 2D numpy.ndarray,
            list(2D numpy.ndarray)}, optional
        The input used to set the flat image. See Calibration module and
        ImageConcersion.convert_image_to_array() for details
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
        need to be tweaked, especially with few co-added images, due to random
        noise. The filtered image is used for the centering algorithms, so for
        a "true test" use kernel_size=1, but be careful, because this may
        lead to needing a fairly high threshold for the noise.
    input_fnum : float
        The focal ratio input into the FCS fiber
    output_fnum : float
        The focal ratio on the output side of the FCS

    Attributes
    ----------
    object_file : string
        File location of the saved object. Only set if self.save_object is
        called
    image_file : string
        File location of the saved image. Only set if self.save_image is called
    data_file : string
        File location of the ImageAnalysis data. Only set if self.save_data is
        called
    image_input : string
        The input for ImageAnalysis instantiation so the calibrated image can
        be created when needed
    calibration : Calibration
        The container for the calibration images
    new_calibration : boolean
        Whether or not self.calibration has been set with new images

    _image_info : ImageInfo
        Container for information about the image (pixel_size, etc.)
    _analysis_info : AnalysisInfo
        Container for information about the analysis
    _frd_info : FRDInfo
        Container for information pertaining to the calculated FRD

    _edges : Edges
        Container for the location of the four fiber edges
    _center : FiberInfo
        Container for the calculated centers of the fiber
    _centroid : FiberInfo
        Container for the calculated centroids of the fiber
    _diameter : FiberInfo
        Container for the calculated diameters of the fiber
    _array_sum : FiberInfo
        Container for the array sums used in the 'circle' methods

    _gaussian_amp : float
        Amplitude of the gaussian fit function
    _gaussian_offset : float
        Offset of the gaussian fit function

    _phi : float
        Golden ratio for the optimization tests
    """
    def __init__(self, image_input, dark=None, ambient=None, flat=None,
                 image_data=None, pixel_size=None, camera=None,
                 magnification=None, threshold=256, kernel_size=9,
                 input_fnum=2.4, output_fnum=2.4):
        # Private attribute initialization
        if image_data is None:
            self._image_info = ImageInfo()
            self._image_info.pixel_size = pixel_size
            self._image_info.camera = camera
            self._image_info.magnification = magnification

            self._analysis_info = AnalysisInfo(kernel_size, threshold)
            self._edges = Edges()
            self._center = FiberInfo('pixel')
            self._centroid = FiberInfo('pixel')
            self._diameter = FiberInfo('value')
            self._array_sum = FiberInfo('value')

            self._frd_info = FRDInfo()
            self._frd_info.input_fnum = input_fnum
            self._frd_info.output_fnum = output_fnum
        else:
            self.load_data(image_data)

        self._phi = (5 ** 0.5 - 1) / 2

        self.object_file = None
        self.image_file = None
        self.data_file = None

        self.image_input = image_input
        self.calibration = Calibration(dark, ambient, flat)
        self.new_calibration = True

        self._gaussian_amp = 0.0
        self._gaussian_offset = 0.0

        self._rectangle_width = 0.0
        self._rectangle_height = 0.0
        self._rectangle_angle = 0.0

        self.set_image_info(image_input)

    #=========================================================================#
    #==== Private Variable Setters ===========================================#
    #=========================================================================#

    def set_image_info(self, image_input):
        """Sets image info using the output of convert_image_to_array()

        Args
        ----
        output_obj : ImageInfo
            The object containing items to add to _image_info

        Sets
        ----
        _image_info.attr
            for all keys in output_obj.__dict__
        """
        _, output_obj = convert_image_to_array(image_input, True)

        for key in output_obj.__dict__:
            if getattr(self._image_info, key) is None:
                setattr(self._image_info, key, getattr(output_obj, key))

        if self._image_info.magnification is None:
            if self._image_info.camera == 'nf' or self._image_info.camera == 'in':
                self._image_info.magnification = 10.0
            else:
                self._image_info.magnification = 1.0

    def set_image_file(self, image_file):
        """Sets the image input.

        Needs to be called if

        Args
        ----
        image_file : string

        Raises
        ------
        RuntimeError
            If the file is not FITS
        """
        if not image_file.endswith('.fit'):
            raise RuntimeError('Please set image file to FITS file')
        self.image_file = image_file

    def set_magnification(self, magnification):
        """Sets the magnification."""
        self._image_info.magnification = magnification

    def set_dark(self, dark):
        """Sets the dark calibration image."""
        self.calibration.dark = dark
        self.new_calibration = True

    def set_ambient(self, ambient):
        """Sets the ambient calibration image."""
        self.calibration.ambient = ambient
        self.new_calibration = True

    def set_flat(self, flat):
        """Sets the flat calibration images."""
        self.calibration.flat = flat
        self.new_calibration = True

    #=========================================================================#
    #==== Saving and Loading Data to File ====================================#
    #=========================================================================#

    def save(self):
        """Save the object, image, and data using the predetermined file names."""
        self.save_object()
        self.save_image()
        self.save_data()

    def save_object(self, file_name=None):
        """Pickle the entire ImageAnalysis object.

        Saves
        -----
        self: ImageAnalysis
            the entire object as .pkl
        """
        if file_name is None and self.object_file is None:
            self.object_file = self._image_info.folder + self.get_camera() + '_object.pkl'
        elif file_name is not None:
            self.object_file = file_name
        save_image_object(self, self.object_file)

    def save_image(self, file_name=None):
        """Save the corrected image as FITS

        Args
        ----
        file_name : {None, string}, optional
            The file name which is used to store the images. The file extension
            should be either '.fit' or '.tif'

        Saves
        -----
        image : 2D numpy.ndarray
            as FITS or TIFF
        """
        if file_name is None and self.image_file is None:
            file_name = self._image_info.folder + self.get_camera() + '_corrected.fit'
        elif file_name is None:
            file_name = self.image_file
        save_array(self.get_image(), file_name)
        self.image_file = file_name

    def save_data(self, file_name=None):
        """Pickle the data and also save the data as a text file dictionary

        Args
        ----
        file_name : {None, string}, optional
            The file name which is used to store the images. The file extension
            should be either '.txt' or '.pkl'

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
        if file_name is None and self.data_file is None:
            file_name = self._image_info.folder + self.get_camera() + '_data.txt'
        elif file_name is None:
            file_name = self.data_file
        save_data(self, file_name)
        self.data_file = file_name

    def load_data(self, file_name=None):
        """Loads data from a text file containing a python dictionary

        Args
        ----
        file_name : string
            The file where the data is located

        Raises
        ------
        RuntimeError
            if the file name does not end in '.txt'

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
        if file_name is None:
            file_name = self.data_file

        if file_name.endswith('.txt'):
            with open(file_name, 'r') as load_file:
                data = literal_eval(load_file.read())
        else:
            raise RuntimeError('Incorrect file type to load into object')

        self._image_info.__dict__ = data['_image_info']
        self._analysis_info.__dict__ = data['_analysis_info']
        self._edges.__dict__ = data['_edges']
        self._center.__dict__ = data['_center']
        self._diameter.__dict__ = data['_diameter']
        self._centroid.__dict__ = data['_centroid']
        self._array_sum.__dict__ = data['_array_sum']

    #=========================================================================#
    #==== Private Variable Getters ===========================================#
    #=========================================================================#

    def get_image(self):
        """Return the corrected image

        This method must be called to get access to the corrected 2D numpy
        array being analyzed. Attempts to access a previously saved image
        under self.image_file or otherwise applies corrections to the raw
        images pulled from their respective files

        Returns
        -------
        image : 2D numpy array
            Image corrected by calibration images
        """
        if self.image_file is not None and not self.new_calibration:
            return convert_image_to_array(self.image_file)
        image = self.calibration.execute_error_corrections(self.get_uncorrected_image(),
                                                           self._image_info)
        self.new_calibration = False
        return image

    def get_uncorrected_image(self):
        """Return the raw image without corrections or filtering.

        Returns
        -------
        uncorrected_image : 2D numpy array
            Raw image or average of images (depending on image_input)
        """
        return convert_image_to_array(self.image_input)

    def get_filtered_image(self, kernel_size=None):
        """Return a median filtered image

        Args
        ----
        kernel_size : {None, int (odd)}, optional
            The side length of the kernel used to median filter the image. Uses
            self._analysis_info.kernel_size if None.

        Returns
        -------
        filtered_image : 2D numpy.ndarray
            The stored image median filtered with the given kernel_size
        None : NoneType
            If image_input is None
        """
        image = self.get_uncorrected_image()
        if image is None:
            return None
        if kernel_size is None:
            kernel_size = self._analysis_info.kernel_size
        filtered_raw_image = filter_image(image, kernel_size)
        return self.calibration.execute_error_corrections(filtered_raw_image,
                                                          self._image_info)

    def get_mesh_grid(self):
        """Return a meshgrid of the same size as the stored image"""
        return mesh_grid_from_array(self.get_image())

    def get_fiber_data(self, method=None, units='microns', **kwargs):
        """Return the fiber center and diameter

        Args
        ----
        method : {None, 'radius', 'gaussian', 'circle', 'edge'}, optional
            The method which is used to calculate the fiber center. If None,
            return the best calculated fiber center in the order 'radius' >
            'gaussian' > 'circle' > 'edge'
        units : {'pixels', 'microns'}, optional
            The units of the returned values
        **kwargs
            The keyworded arguments to pass to the centering method

        Sets
        ----
        _diameter.method : float
            The diameter of the fiber face in the context of the given method
        _center.method : {'x': float, 'y': float}
            The center of the fiber face in the context of the given method

        Returns
        -------
        _center.method.y : float
            in the given units
        _center.method.x : float
            in the given units
        _diameter.method : float
            in the given units
        """
        center_y, center_x = self.get_fiber_center(method, units=units, **kwargs)
        kwargs['show_image'] = False
        diameter = self.get_fiber_diameter(method, units=units, **kwargs)
        return center_y, center_x, diameter

    def get_fiber_radius(self, method=None, units='pixels', **kwargs):
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
        **kwargs
            The keyworded arguments to pass to the centering method

        Sets
        ----
        _diameter.method : float
            The diameter of the fiber face in the context of the given method
        _center.method : {'x': float, 'y': float}
            The center of the fiber face in the context of the given method

        Returns
        -------
        _diameter.method / 2.0 : float
            in the given units
        """
        return self.get_fiber_diameter(method, units=units, **kwargs) / 2.0

    def get_fiber_diameter(self, method=None, units='pixels', **kwargs):
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
        **kwargs
            The keyworded arguments to pass to the centering method

        Sets
        ----
        _diameter.method : float
            The diameter of the fiber face in the context of the given method
        _center.method : {'x': float, 'y': float}
            The center of the fiber face in the context of the given method

        Returns
        -------
        _diameter.method : float
            in the given units
        """
        if method is None:
            if self._image_info.camera != 'ff':
                if self._diameter.radius is not None:
                    method = 'radius'
                elif self._diameter.circle is not None:
                    method = 'circle'
                else:
                    method = 'edge'
            elif self._diameter.gaussian is not None:
                method = 'gaussian'
            else:
                method = 'edge'

        if getattr(self._diameter, method) is None:
            self.set_fiber_diameter(method, **kwargs)

        diameter = getattr(self._diameter, method)

        return self.convert_pixels_to_units(diameter, units)

    def get_fiber_center(self, method=None, units='pixels', **kwargs):
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
        **kwargs
            The keyworded arguments to pass to the centering method

        Sets
        ----
        _diameter.method : float
            The diameter of the fiber face in the context of the given method
        _center.method : {'x': float, 'y': float}
            The center of the fiber face in the context of the given method

        Returns
        -------
        _center.method.y : float
            in the given units
        _center.method.x : float
            in the given units
        """
        if method is None:
            if self._image_info.camera != 'ff':
                if self._center.radius.x is not None:
                    method = 'radius'
                elif self._center.circle.x is not None:
                    method = 'circle'
                else:
                    method = 'edge'
            elif self._center.gaussian.x is not None:
                method = 'gaussian'
            else:
                method = 'edge'

        if getattr(self._center, method).x is None or method == 'circle':
            self.set_fiber_center(method, **kwargs)

        center = (getattr(self._center, method).y,
                  getattr(self._center, method).x)

        return self.convert_pixels_to_units(center, units)

    def get_fiber_centroid(self, method=None, units='pixels', **kwargs):
        """Getter for the fiber centroid

        Args
        ----
        method : {None, 'full', 'edge', 'radius', 'gaussian', 'circle'}, optional
            See set_fiber_centroid() for method details. If no method is given,
            chooses the most precise method already calculated in the order
            'radius' > 'gaussian' > 'circle' > 'edge' > 'full'
        units : {'pixels', 'microns'}, optional
            The units of the returned values
        **kwargs
            The keyworded arguments to pass to the centering and centroiding
            methods

        Sets
        ----
        _centroid.method : {'x': float, 'y': float}
            The centroid of the fiber face in the context of the given method

        Returns
        -------
        _centroid.method.y : float
            in the given units
        _centroid.method.x : float
            in the given units
        """
        if method is None:
            if self._image_info.camera != 'ff':
                if self._centroid.radius.x is not None:
                    method = 'radius'
                elif self._centroid.circle.x is not None:
                    method = 'circle'
                else:
                    method = 'edge'
            elif self._centroid.gaussian.x is not None:
                method = 'gaussian'
            elif self._centroid.edge.x is not None:
                method = 'edge'
            else:
                method = 'full'

        if getattr(self._centroid, method).x is None:
            self.set_fiber_centroid(method, **kwargs)
        centroid = (getattr(self._centroid, method).y,
                    getattr(self._centroid, method).x)

        return self.convert_pixels_to_units(centroid, units)

    def get_rectangle_fit(self):
        """Return the best rectangle fit for the image"""
        if self._center.rectangle.x is None:
            self.set_fiber_center(method='rectangle')
        rectangle_fit = rectangle_array(self.get_mesh_grid(),
                                        self._center.rectangle.x,
                                        self._center.rectangle.y,
                                        self._rectangle_width,
                                        self._rectangle_height,
                                        self._rectangle_angle
                                       ).reshape(self.get_height(),
                                                 self.get_width())
        return rectangle_fit

    def get_gaussian_fit(self):
        """Return the best gaussian fit for the image

        Returns
        -------
        _fit.gaussian : 2D numpy.ndarray
        """
        if self._center.gaussian.x is None:
            self.set_fiber_center(method='gaussian')
        gaussian_fit = gaussian_array(self.get_mesh_grid(),
                                      self._center.gaussian.x,
                                      self._center.gaussian.y,
                                      self._diameter.gaussian / 2.0,
                                      self._gaussian_amp,
                                      self._gaussian_offset
                                     ).reshape(self.get_height(),
                                               self.get_width())

        if self._image_info.camera == 'in':
            y0, x0 = self.get_fiber_center(method='edge')
            radius = self.get_fiber_radius(method='edge')
            filtered_image = self.get_filtered_image()
            fiber_face = circle_array(self.get_mesh_grid(), x0, y0, radius)
            fiber_face *= np.median(crop_image(filtered_image, x0, y0, radius)[0])
            gaussian_fit += fiber_face

        return gaussian_fit

    def get_polynomial_fit(self, deg=6, x0=None, y0=None):
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
        if y0 is None or x0 is None:
            y0, x0 = self.get_fiber_center()
        return polynomial_fit(self.get_image(), deg, x0, y0)

    def get_tophat_fit(self):
        """Return the circle array that best covers the fiber face

        Returns
        -------
        circle_array : numpy.ndarray (2D)
            Circle array centered at best calculated center and with best
            calculated diameter
        """
        y0, x0 = self.get_fiber_center()
        radius = self.get_fiber_radius()
        return self.get_image().max() * circle_array(self.get_mesh_grid(), x0, y0, radius, res=1)

    def get_image_info(self, info_type=None):
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
        elif info_type in self._image_info.__dict__:
            if getattr(self._image_info, info_type) is None:
                raise RuntimeError(info_type + ' needs to be set externally')
            return getattr(self._image_info, info_type)
        raise RuntimeError('Incorrect string for image info property')

    def get_analysis_info(self, info_type=None):
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
        elif info_type in self._analysis_info.__dict__:
            if getattr(self._analysis_info, info_type) is None:
                raise RuntimeError(info_type + ' needs to be set externally')
            return getattr(self._analysis_info, info_type)
        raise RuntimeError('Incorrect string for image info property')

    def get_height(self, units='pixels'):
        """Return the image height in units"""
        return self.convert_pixels_to_units(self._image_info.height, units)

    def get_width(self, units='pixels'):
        """Return the image width in units"""
        return self.convert_pixels_to_units(self._image_info.width, units)

    def get_magnification(self):
        """Return the magnification"""
        if self._image_info.magnification is None:
            raise RuntimeError('Magnification needs to be set externally')
        return self._image_info.magnification

    def get_camera(self):
        """Return the string denoting the camera type"""
        if self._image_info.camera is None:
            raise RuntimeError('Camera needs to be set externally')
        return self._image_info.camera

    def get_pixel_size(self):
        """Return the pixel size in microns"""
        if self._image_info.pixel_size is None:
            raise RuntimeError('Pixel Size needs to be set externally')
        return self._image_info.pixel_size

    def get_num_images(self):
        """Return the number of images that created the object"""
        if self._image_info.num_images is None:
            raise RuntimeError('Number of images improperly set')
        return self._image_info.num_images

    def get_input_fnum(self):
        """Return the focal ratio of the FCS input side."""
        return self._frd_info.input_fnum

    def get_output_fnum(self):
        """Return the focal ratio of the FCS output side."""
        return self._frd_info.output_fnum

    def get_frd_info(self, new=False, **kwargs):
        """Return the FRDInfo object and sets it where appropriate.

        See set_frd_info for more information on kwargs

        Args
        ----
        new : boolean
            If new is True, recalculate the FRD info

        Returns
        -------
        self._frd_info : FRDInfo
            Container for FRD information. See Containers.FRDInfo
        """
        if new or not self._frd_info.energy_loss:
            self.set_frd_info(**kwargs)
        return self._frd_info

    #=========================================================================#
    #==== Calculating Output F Number ========================================#
    #=========================================================================#

    def set_frd_info(self, f_lim=(2.3, 6.0), res=0.1, fnum_diameter=0.95):
        """Calculate the encircled energy for various focal ratios

        Args
        ----
        f_lim : (float, float)
            the limits of the focal ratio calculations
        res : float
            the spacing between each focal ratio when calculating encircled
            energy
        fnum_diameter : float
            the fraction of the total encircled energy at which the output
            focal ratio is set

        Sets
        ----
        _frd_info.output_fnum : float
            the approximate focal ratio inside which 95% of the total encircled energy
            is included
        _frd_info.energy_loss : float
            the loss of energy when the output focal ratio equals the input focal ratio
            given as a percent
        _frd_info.encircled_energy_fnum : list(float)
            list of the focal ratios used to calculate encircled energy
        _frd_info.encircled_energy : list(float)
            list of the encircled energy at each given focal ratio
        """
        center_y, center_x = self.get_fiber_centroid(method='full')

        fnums = list(np.arange(f_lim[0], f_lim[1] + res, res))
        energy_loss = None
        output_fnum = None
        encircled_energy = []
        image = self.get_image()
        for fnum in fnums:
            radius = self.convert_fnum_to_radius(fnum, units='pixels')
            isolated_circle = isolate_circle(image,
                                             center_x,
                                             center_y,
                                             radius)
            iso_circ_sum = sum_array(isolated_circle)
            encircled_energy.append(iso_circ_sum)
            if abs(fnum - self._frd_info.input_fnum) < res / 2.0:
                energy_loss = 100 * (1 - iso_circ_sum / encircled_energy[0])
            if iso_circ_sum / encircled_energy[0] >= fnum_diameter:
                output_fnum = fnum

        self._frd_info.output_fnum = output_fnum
        self._frd_info.energy_loss = energy_loss
        self._frd_info.encircled_energy_fnum = fnums
        self._frd_info.encircled_energy = list(np.array(encircled_energy) / encircled_energy[0])

    #=========================================================================#
    #==== Image Centroiding ==================================================#
    #=========================================================================#

    def set_fiber_centroid(self, method='full', radius_factor=1.0,
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
        _centroid.method : Pixel
            The centroid of the image in the context of the given method
        """
        image = self.get_filtered_image()
        if method == 'full':
            image_array_iso = image
        else:
            y0, x0 = self.get_fiber_center(method=method,
                                           show_image=False,
                                           **kwargs)
            radius = self.get_fiber_radius(method=method,
                                           show_image=False,
                                           **kwargs)
            image_array_iso = isolate_circle(image, x0, y0,
                                             radius*radius_factor, res=1)

        x_array, y_array = self.get_mesh_grid()
        getattr(self._centroid, method).x = ((image_array_iso * x_array).sum()
                                             / image_array_iso.sum())
        getattr(self._centroid, method).y = ((image_array_iso * y_array).sum()
                                             / image_array_iso.sum())

        if show_image:
            if method == 'gaussian':
                plot_overlaid_cross_sections(image, self.get_gaussian_fit(),
                                             self._center.gaussian.y,
                                             self._center.gaussian.x)
            plot_dot(image,
                     getattr(self._centroid, method).y,
                     getattr(self._centroid, method).x)
            show_plots()

    #=========================================================================#
    #==== Image Centering ====================================================#
    #=========================================================================#

    def set_fiber_data(self, method, **kwargs):
        """Set the fiber center, diameter, and centroid using the same method

        Args
        ----
        method : {'edge', 'radius', 'gaussian', 'circle'}
            Uses the respective method to find the fiber center
        **kwargs
            The keyworded arguments to pass to the centering method

        Sets
        ----
        _centroid.method : Pixel
            The centroid of the image in the context of the given method
        _center.method : Pixel
            The center of the fiber face in the context of the given method
        _diameter.method : float
            The diameter of the fiber face in the context of the given method
        """
        self.set_fiber_center(method, **kwargs)
        self.set_fiber_centroid(method, **kwargs)

    def set_fiber_diameter(self, method, **kwargs):
        """Set the fiber diameter using given method

        Args
        ----
        method : {'edge', 'radius', 'gaussian', 'circle'}
            Uses the respective method to find the fiber center
        **kwargs :
            The keyworded arguments to pass to the centering method

        Sets
        ----
        _diameter.method : float
            The diameter of the fiber face in the context of the given method
        _center.method : Pixel
            The center of the fiber face in the context of the given method

        Raises
        ------
        RuntimeError
            cannot accept the 'circle' method when setting the diameter since
            it requires a known radius to run
        """
        if method == 'circle':
            raise RuntimeError('Fiber diameter cannot be set by circle method')
        self.set_fiber_center(method, **kwargs)

    def set_fiber_center(self, method, show_image=False, **kwargs):
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
        if method == 'radius':
            self.set_fiber_center_radius_method(**kwargs)
        elif method == 'edge':
            self.set_fiber_center_edge_method()
        elif method == 'circle':
            self.set_fiber_center_circle_method(**kwargs)
        elif method == 'gaussian':
            self.set_fiber_center_gaussian_method()
        elif method == 'rectangle':
            self.set_fiber_center_rectangle_method(**kwargs)
        else:
            raise RuntimeError('Incorrect string for fiber centering method')

        if show_image:
            x0 = getattr(self._center, method).x
            y0 = getattr(self._center, method).y
            r = getattr(self._diameter, method) / 2.0
            image = self.get_filtered_image()

            if method == 'gaussian':
                plot_overlaid_cross_sections(image, self.get_gaussian_fit(),
                                             y0, x0)
                plot_dot(image, y0, x0)
                show_plots()
            else:
                plot_image_array(remove_circle(image, x0, y0, r, res=1))
                plot_overlaid_cross_sections(image, image.max() / 2.0
                                             * circle_array(self.get_mesh_grid(),
                                                            x0, y0, r, res=1),
                                             y0, x0)
                show_plots()

    def set_fiber_center_rectangle_method(self, radius=None, **kwargs):
        """Set fiber center using a rectangle mask

        Uses Scipy.optimize.curve_fit method to fit fiber image to 
        rectangle_array().

        Sets
        ----
        _diameter.rectangle : float
            sqrt(height^2 + width^2) of the rectangle)
        _center.rectangle : Pixel()
            Center of the fiber in the rectangle method context
        _fit.rectangle : 2D numpy.ndarray
            Best rectangle fit for the fiber image
        """
        if self._edges.left is None:
            self.set_fiber_edges()

        image = self.get_filtered_image()

        left = np.array([image[:, self._edges.left].argmax(),
                         self._edges.left])
        right = np.array([image[:, self._edges.right].argmax(),
                         self._edges.right])
        top = np.array([self._edges.top,
                        image[self._edges.top, :].argmax()])
        bottom = np.array([self._edges.bottom,
                           image[self._edges.bottom, :].argmax()])

        radius = (np.sqrt(((right - left)**2).sum())
                  + np.sqrt(((bottom - top)**2).sum())) / 4.0

        self.set_fiber_center_circle_method(radius, **kwargs)

        self._center.rectangle.x = self._center.circle.x
        self._center.rectangle.y = self._center.circle.y
        self._diameter.rectangle = radius * 2.0

        # fiber_y0, fiber_x0 = self.get_fiber_center(method='edge')
        # fiber_width = np.sqrt(self.get_fiber_diameter(method='edge')**2 / 2.0)
        # filtered_image = self.get_filtered_image()

        # initial_guess = (fiber_x0, fiber_y0, fiber_width, fiber_width, 0.0)
        # _, opt_parameters = rectangle_fit(filtered_image,
        #                                   initial_guess=initial_guess,
        #                                   full_output=True)
        # print opt_parameters
        # self._center.rectangle.x = opt_parameters[0]
        # self._center.rectangle.y = opt_parameters[1]
        # self._rectangle_width = opt_parameters[2]
        # self._rectangle_height = opt_parameters[3]
        # self._diameter.rectangle = np.sqrt(self._rectangle_width**2 + self._rectangle_height**2)
        # self._rectangle_angle = opt_parameters[4]

    def set_fiber_center_gaussian_method(self):
        """Set fiber center using a Gaussian Fit

        Uses Scipy.optimize.curve_fit method to fit fiber image to
        gaussian_array(). The radius found extends to 2-sigma of the gaussian
        therefore encompassing ~95% of the imaged light. Use previous methods
        of center-finding to approximate the location of the center

        Sets
        ----
        _diameter.gaussian : float
            Diameter of the fiber in the gaussian method context
        _center.gaussian : {'x': float, 'y': float}
            Center of the fiber in the gaussian method context
        _fit.gaussian : 2D numpy.ndarray
            Best gaussian fit for the fiber image
        """
        fiber_y0, fiber_x0 = self.get_fiber_center(method='edge')
        fiber_radius = self.get_fiber_radius(method='edge')
        filtered_image = self.get_filtered_image()

        if self._image_info.camera == 'in':
            radius = -1
            factor = 0.9
            fiber_face = circle_array(self.get_mesh_grid(), fiber_x0, fiber_y0,
                                      fiber_radius, res=1)
            fiber_face *= np.median(crop_image(filtered_image, fiber_x0,
                                               fiber_y0, fiber_radius)[0])
            while radius < 1 or radius > fiber_radius:
                factor += 0.1
                cropped_image, new_x0, new_y0 = crop_image(filtered_image,
                                                           fiber_x0, fiber_y0,
                                                           fiber_radius*factor)

                initial_guess = (new_x0, new_y0, 100 / self.get_pixel_size(),
                                 cropped_image.max(), cropped_image.min())
                try:
                    fit, opt_parameters = gaussian_fit(cropped_image,
                                                       initial_guess=initial_guess,
                                                       full_output=True)

                    radius = abs(opt_parameters[2])
                except RuntimeError:
                    radius = -1

            x0 = opt_parameters[0] + int(fiber_x0-fiber_radius*factor)
            y0 = opt_parameters[1] + int(fiber_y0-fiber_radius*factor)
            amp = opt_parameters[3]
            offset = opt_parameters[4]

        else:
            initial_guess = (fiber_x0, fiber_y0, fiber_radius,
                             filtered_image.max(), filtered_image.min())

            _, opt_parameters = gaussian_fit(filtered_image,
                                             initial_guess=initial_guess,
                                             full_output=True)
            x0 = opt_parameters[0]
            y0 = opt_parameters[1]
            radius = abs(opt_parameters[2])
            amp = opt_parameters[3]
            offset = opt_parameters[4]

        self._center.gaussian.x = x0
        self._center.gaussian.y = y0
        self._diameter.gaussian = radius * 2.0
        self._gaussian_amp = amp
        self._gaussian_offset = offset

    def set_fiber_center_radius_method(self, tol=.03, test_range=None):
        """Set fiber center using dark circle with varying radius

        Uses a golden mean optimization method to find the optimal radius of the
        dark circle that covers the fiber image used in
        get_fiber_centerCircleMethod(). The optimization is for a parameter
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
        _diameter.radius : float
            Diameter of the fiber in the radius method context
        _center.radius : {'x': float, 'y': float}
            Center of the fiber in the radius method context
        _diameter.circle : float
            Also uses the circle method, therefore changes this value
        _center.circle : float
            Also uses the circle method, therefore chnages this value
        """        
        image = self.get_filtered_image()

        # Initialize range of tested radii
        r = np.zeros(4).astype(float)

        if test_range is not None:
            approx_radius = self.get_fiber_radius()
            test_range = test_range / 2.0

            r[0] = approx_radius - test_range
            if r[0] < 0.0:
                r[0] = 0.0
            r[3] = approx_radius + test_range
        else:
            r[0] = 0
            r[3] = min(self._image_info.height, self._image_info.width) / 2.0

        r[1] = r[0] + (1 - self._phi) * (r[3] - r[0])
        r[2] = r[0] + self._phi * (r[3] - r[0])

        array_sum = np.zeros(2).astype(float)
        for i in xrange(2):
            self.set_fiber_center(method='circle', radius=r[i+1],
                                  tol=tol, test_range=test_range, image=image)
            array_sum[i] = (self._array_sum.circle
                            + self._analysis_info.threshold
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

            self.set_fiber_center(method='circle', radius=r[min_index+1],
                                  tol=tol, test_range=test_range, image=image)
            array_sum[min_index] = (self._array_sum.circle
                                    + self._analysis_info.threshold
                                    * np.pi * r[min_index+1]**2)

            min_index = np.argmin(array_sum) # Integer 0 or 1 for min of r[1], r[2]

        self._diameter.radius = r[min_index+1] * 2
        self._center.radius.y = self._center.circle.y
        self._center.radius.x = self._center.circle.x
        self._array_sum.radius = np.amin(array_sum)

    def set_fiber_center_circle_method(self, radius, tol=.03, test_range=None, image=None):
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
        image : 2d numpy.ndarray, optional
            The image being analyzed. This is only useful for the radius_method.
            Probably not for use outside the class.

        Sets
        ----
        _diameter.circle : float
            Diameter of the fiber in the circle method context
        _center.circle : {'x': float, 'y': float}
            Center of the fiber in the circle method context
        _diameter.edge : float
            If test_range is not None, approximates the circle's center using
            the edge method
        _center.edge : float
            If test_range is not None, approximates the circle's center using
            the edge method
        """
        res = int(1.0/tol)
        if image is None:
            image = self.get_filtered_image()

        # Create four "corners" to test center of the removed circle
        x = np.zeros(4).astype(float)
        y = np.zeros(4).astype(float)

        if test_range is not None:
            approx_center = self.get_fiber_center(method='edge', show_image=False)
            test_range = test_range / 2.0

            x[0] = approx_center[1] - test_range
            if x[0] < radius:
                x[0] = radius
            x[3] = approx_center[1] + test_range
            if x[3] > self._image_info.width - radius:
                x[3] = self._image_info.width - radius

            y[0] = approx_center[0] - test_range
            if y[0] < radius:
                y[0] = radius
            y[3] = approx_center[0] + test_range
            if y[3] > self._image_info.height - radius:
                y[3] = self._image_info.height - radius

        else:
            x[0] = radius
            x[3] = self._image_info.width - radius

            y[0] = radius
            y[3] = self._image_info.height - radius

        x[1] = x[0] + (1 - self._phi) * (x[3] - x[0])
        x[2] = x[0] + self._phi * (x[3] - x[0])

        y[1] = y[0] + (1 - self._phi) * (y[3] - y[0])
        y[2] = y[0] + self._phi * (y[3] - y[0])

        # Initialize array sums to each corner
        array_sum = np.zeros((2, 2)).astype(float)
        for i in xrange(2):
            for j in xrange(2):
                removed_circle_array = remove_circle(image,
                                                     x[i+1], y[j+1],
                                                     radius, res=1)
                array_sum[j, i] = sum_array(removed_circle_array)

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
                        temp_res = res
                        if abs(x[3] - x[0]) < 10*tol or abs(y[3] - y[0]) < 10*tol:
                            temp_res = 1
                        removed_circle_array = remove_circle(image,
                                                             x[i+1], y[j+1],
                                                             radius, temp_res)
                        array_sum[j, i] = sum_array(removed_circle_array)

            min_index = np.unravel_index(np.argmin(array_sum), (2, 2))

        self._center.circle.x = x[min_index[1]+1]
        self._center.circle.y = y[min_index[0]+1]
        self._diameter.circle = radius * 2.0
        self._array_sum.circle = np.amin(array_sum)

    def set_fiber_center_edge_method(self):
        """TAverages the fiber edges to set the fiber center

        Sets
        ----
        self._center.edge.y : float
        self._center.edge.x : float
        """
        self.set_fiber_edges()

        self._center.edge.y = (self._edges.top + self._edges.bottom) / 2.0
        self._center.edge.x = (self._edges.left + self._edges.right) / 2.0

    def set_fiber_edges(self):
        """Set fiber edge pixel values

        Sets the left, right, top, and bottom edges of the fiber by finding where
        the maxima of each row and column cross the given threshold. Also sets
        the width of the fiber by the maximum of the horizontal and vertical
        lengths

        Sets
        ----
        self._edges.left : float
        self._edges.right : float
        self._edges.top : float
        self._edges.bottom : float
        self._diameter.edge : float
        """
        filtered_image = self.get_filtered_image()

        left = -1
        right = -1
        for index in xrange(self._image_info.width):
            if left < 0:
                if filtered_image[:, index].max() > self._analysis_info.threshold:
                    left = index
            else:
                if filtered_image[:, index].max() > self._analysis_info.threshold:
                    right = index

        top = -1
        bottom = -1
        for index in xrange(self._image_info.height):
            if top < 0:
                if filtered_image[index, :].max() > self._analysis_info.threshold:
                    top = index
            else:
                if filtered_image[index, :].max() > self._analysis_info.threshold:
                    bottom = index

        self._edges.left = left
        self._edges.right = right
        self._edges.top = top
        self._edges.bottom = bottom
        self._diameter.edge = ((right - left) + (bottom - top)) / 2.0

    #=========================================================================#
    #==== Overriding Methods =================================================#
    #=========================================================================#

    def convert_pixels_to_units(self, value, units):
        """Returns the value in the proper units"""
        return convert_pixels_to_units(value,
                                       self.get_pixel_size(),
                                       self.get_magnification(),
                                       units)

    def convert_fnum_to_radius(self, fnum, units):
        """Returns the value in the proper units"""
        return convert_fnum_to_radius(fnum,
                                      self.get_pixel_size(),
                                      self.get_magnification(),
                                      units)

    def plot_cross_sections(self, image_array=None, row=None, column=None):
        """Plots cross sections across the center of the fiber"""
        if image_array is None:
            image_array = self.get_image()
        if row is None:
            row = self.get_fiber_center()[0]
        if column is None:
            column = self.get_fiber_center()[1]
        plot_cross_sections(image_array, row, column)

    def show_image_array(self, image_array=None):
        """Shows the calibrated image"""
        if image_array is None:
            image_array = self.get_image()
        show_image_array(image_array)


if __name__ == "__main__":
    from fiber_properties.input_output import load_image_object

    folder = 'C:/Libraries/Box Sync/ExoLab/Fiber_Characterization/Image Analysis/data/scrambling/2016-08-05 Prototype Core Extension 1/'

    images = [folder + 'Shift_00/in_' + str(i).zfill(3) + '.fit' for i in xrange(10)]
    dark = [folder + 'Dark/in_' + str(i).zfill(3) + '.fit' for i in xrange(10)]
    ambient = [folder + 'Ambient/in_' + str(i).zfill(3) + '.fit' for i in xrange(10)]

    im_obj = ImageAnalysis(images, dark, ambient)

    tol = 0.25
    test_range = 5
    factor = 1.0

    im_obj.show_image_array()
    for key in im_obj.get_image_info().__dict__:
        print key + ': ' + str(im_obj.get_image_info().__dict__[key])
    for key in im_obj.get_analysis_info().__dict__:
        print key + ': ' + str(im_obj.get_analysis_info().__dict__[key])
    print
    print 'Centroid:'
    print im_obj.get_fiber_centroid(method='full', radius_factor=factor)
    print
    print 'Edge:'
    print 'center:', im_obj.get_fiber_center(method='edge', show_image=False)
    print 'centroid:', im_obj.get_fiber_centroid(method='edge', radius_factor=factor)
    print 'diameter:', im_obj.get_fiber_diameter(method='edge', units='microns'), 'microns'
    print
    print 'Radius:'
    print 'center:', im_obj.get_fiber_center(method='radius', tol=tol, test_range=test_range, show_image=False)
    print 'centroid:', im_obj.get_fiber_centroid(method='radius', radius_factor=factor)
    print 'diameter:', im_obj.get_fiber_diameter(method='radius', units='microns'), 'microns'
    print
    print 'Gaussian:'
    print 'center:', im_obj.get_fiber_center(method='gaussian', show_image=False)
    print 'centroid:', im_obj.get_fiber_centroid(method='gaussian', radius_factor=factor)
    print 'diameter:', im_obj.get_fiber_diameter(method='gaussian', units='microns'), 'microns'

    im_obj.save()

    new_im_obj = load_image_object(im_obj.object_file, im_obj.image_file)

    new_im_obj.show_image_array()
    for key in new_im_obj.get_image_info().__dict__:
        print key + ': ' + str(im_obj.get_image_info().__dict__[key])
    for key in new_im_obj.get_analysis_info().__dict__:
        print key + ': ' + str(im_obj.get_analysis_info().__dict__[key])
    print
    print 'Centroid:'
    print new_im_obj.get_fiber_centroid(method='full', radius_factor=factor)
    print
    print 'Edge:'
    print 'center:', new_im_obj.get_fiber_center(method='edge')
    print 'centroid:', new_im_obj.get_fiber_centroid(method='edge', radius_factor=factor)
    print 'diameter:', im_obj.get_fiber_diameter(method='edge', units='microns'), 'microns'
    print
    print 'Radius:'
    print 'center:', new_im_obj.get_fiber_center(method='radius', tol=tol, test_range=test_range)
    print 'centroid:', new_im_obj.get_fiber_centroid(method='radius', radius_factor=factor)
    print 'diameter:', im_obj.get_fiber_diameter(method='radius', units='microns'), 'microns'
    print
    print 'Gaussian:'
    print 'center:', new_im_obj.get_fiber_center(method='gaussian')
    print 'centroid:', new_im_obj.get_fiber_centroid(method='gaussian', radius_factor=factor)
    print 'diameter:', im_obj.get_fiber_diameter(method='gaussian', units='microns'), 'microns'
