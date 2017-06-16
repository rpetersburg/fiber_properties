"""containers.py was written by Ryan Petersburg for use with fiber
characterization on the EXtreme PREcision Spectrograph

The classes in this module are used as containers for information (similar to
dictionaries) in the FiberImage class and image_conversion.py functions.
These are used instead of dictionaries due to the simplicity of attribute
instantiation so that the information is ALWAYS either a value or NONE rather
than an empty slot in a dictionary.
"""
from collections import Iterable
import numpy as np

#=============================================================================#
#===== Metadata Containers ===================================================#
#=============================================================================#

class Edges(object):
    """Container for the fiber image edges."""
    def __init__(self, **kwargs):
        self.left = Pixel(**kwargs)
        self.right = Pixel(**kwargs)
        self.top = Pixel(**kwargs)
        self.bottom = Pixel(**kwargs)

    def __iter__(self):
        for corner in [self.left, self.top, self.right, self.bottom]:
            yield corner

class FiberInfo(object):
    """Container for information concerning the fiber grouped by method."""
    def __init__(self, info=None, **kwargs):
        if info == 'pixel':
            self.edge = Pixel(**kwargs)
            self.radius = Pixel(**kwargs)
            self.circle = Pixel(**kwargs)
            self.gaussian = Pixel(**kwargs)
            self.rectangle = Pixel(**kwargs)
            self.full = Pixel()
        elif info == 'value':
            self.edge = None
            self.radius = None
            self.circle = None
            self.gaussian = None
            self.rectangle = None
            self.full = None

class RectangleInfo(object):
    """Container for information about a rectangle"""
    def __init__(self):
        self.height = None
        self.width = None
        self.angle = None

class Pixel(object):
    """Container for the x and y position of a pixel."""
    def __init__(self, x=None, y=None, units='pixels',
                 pixel_size=None, magnification=None):
        self._x = x
        self._y = y
        self.units = units
        self.magnification = magnification
        self.pixel_size = pixel_size

    def __repr__(self):
        return '(' + str(self.y) + ', ' + str(self.x) + ')'

    def __add__(self, pixel):
        self.test_equality(pixel)
        return Pixel(self.x + pixel.x, self.y + pixel.y, self.units,
                     self.pixel_size, self.magnification)

    def __sub__(self, pixel):
        self.test_equality(pixel)
        return Pixel(self.x - pixel.x, self.y - pixel.y, self.units,
                     self.pixel_size, self.magnification)

    def convert_pixels_to_units(self, value, units):
        return convert_pixels_to_units(value, self.pixel_size,
                                       self.magnification, units)

    def as_tuple(self):
        return (self._x, self._y)

    def as_array(self):
        return np.array(self.as_tuple())

    def test_equality(self, pixel):
        if self.units != pixel.units:
            raise RuntimeError('Two Pixel objects must have same units')
        if self.magnification != pixel.magnification:
            raise RuntimeError('Two Pixel objects must have same magnification')

    def set_pixel(self, *args):
        if len(args) == 1:
            self.x = args[0][0]
            self.y = args[0][1]
        elif len(args) == 2:
            self.x = args[0]
            self.y = args[1]
        else:
            raise RuntimeError('Pixel can only be set by tuple or two values')

    @property
    def x(self):
        return self.convert_pixels_to_units(self._x, self.units)

    @x.setter
    def x(self, x):
        if self.units != 'pixels':
            raise RuntimeError('Please do not set pixel value in microns')
        self._x = x

    @property
    def y(self):
        return self.convert_pixels_to_units(self._y, self.units)

    @y.setter
    def y(self, y):
        if self.units != 'pixels':
            raise RuntimeError('Please do not set pixel value in microns')
        self._y = y

#=============================================================================#
#===== Useful Functions ======================================================#
#=============================================================================#

def convert_pixels_to_units(value, pixel_size, magnification, units):
    """Converts a value or iterable from pixels to given units"""
    if units == 'pixels':
        return value
    elif units == 'microns':
        if isinstance(value, Pixel):
            pixel = Pixel(value.x, value.y)
            pixel.pixel_size = pixel_size
            pixel.magnification = magnification
            pixel.units = units
            return pixel
        elif isinstance(value, Iterable):
            return tuple(np.array(value) * pixel_size / magnification)
        return value * pixel_size / magnification
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
    else:
        raise RuntimeError('Incorrect string for units')

#=============================================================================#
#===== Fiber Property Containers =============================================#
#=============================================================================#

class FRDInfo(object):
    """Container for FRD information

    Attributes
    ----------
    input_fnum : float or list(float)
        list of the given input focal ratios
    encircled_energy : list(float) or list(list(float))
        list of encircled energies for each input_fnum
    encircled_energy_fnum : list(float) or list(list(float))
        independent variable (output f/#) corresponding to each
        encircled energy
    energy_loss : float or list(float)
        energy loss for each input focal ratio
    output_fnum : float or list(float)
        calculated output focal ratio for each input focal ratio
    """
    def __init__(self):
        self.input_fnum = []
        self.output_fnum = []
        self.encircled_energy_fnum = []
        self.encircled_energy = []
        self.energy_loss = []

class ScramblingInfo(object):
    """Container for scrambling gain information

    Attributes
    ----------
    in_x : list(float)
        List of the input centroid x positions
    in_y : list(float)
        List of the input centroid y positions
    out_x : list(float)
        List of the output centroid x positions
    out_y : list(float)
        List of the output centroid y positions
    scrambling_gain : list(float)
        List of the calculated scrambling gains
    in_d : list(float)
        List of all possible permutations of input shift distances
    out_d : list(float)
        List of the resultant output centroid shifts due to in_d
    """
    def __init__(self):
        self.in_x = []
        self.in_y = []
        self.out_x = []
        self.out_y = []
        self.scrambling_gain = []
        self.in_d = []
        self.out_d = []

class ModalNoiseInfo(object):
    """Container for modal noise information

    Attributes
    ----------
    tophat : float
    polynomial: float
    gaussian: float
    gradient : float
    contrast : float
    gini : float
    entropy : float
    fft : FFTInfo
    """
    def __init__(self):
        self.tophat = None
        self.polynomial = None
        self.gaussian = None
        self.gradient = None
        self.contrast = None
        self.gini = None
        self.entropy = None
        self.fft = None

class FFTInfo(object):
    """Container for fast fourier transform information

    Attributes
    ----------
    power : list(float)
        normalized, azimuthally averaged power spectrum
    freq : list(float)
        respective frequencies in 1/um
    """
    def __init__(self, power=None, freq=None):
        self.power = power
        self.freq = freq