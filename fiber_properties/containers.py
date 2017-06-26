"""containers.py was written by Ryan Petersburg for use with fiber
characterization on the EXtreme PREcision Spectrograph

The classes in this module are used as containers for information (similar to
dictionaries) in the FiberImage class and image_conversion.py functions.
These are used instead of dictionaries due to the simplicity of attribute
instantiation so that the information is ALWAYS either a value or NONE rather
than an empty slot in a dictionary.
"""
from __future__ import division
from collections import Iterable
import numpy as np

#=============================================================================#
#===== Metadata Containers ===================================================#
#=============================================================================#

class Edges(object):
    """Container for the fiber image edges."""
    def __init__(self):
        self.left = Pixel()
        self.right = Pixel()
        self.top = Pixel()
        self.bottom = Pixel()

    def __iter__(self):
        for corner in [self.left, self.top, self.right, self.bottom]:
            yield corner

class FiberInfo(object):
    """Container for information concerning the fiber grouped by method."""
    def __init__(self, info=None):
        if info == 'pixel':
            self.edge = Pixel()
            self.radius = Pixel()
            self.circle = Pixel()
            self.gaussian = Pixel()
            self.rectangle = Pixel()
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
    def __init__(self, x=None, y=None, units='pixels'):
        self._x = x
        self._y = y
        self.units = units

    def __repr__(self):
        return '(' + str(self.x) + ', ' + str(self.y) + ')'

    def __add__(self, other):
        if isinstance(other, Pixel):
            x = self._x + other.x
            y = self._y + other.y
        elif isinstance(other, Iterable):
            return self + Pixel(*other)
        else:
            x = self._x + other
            y = self._y + other
        return Pixel(x,y)

    __radd__ = __add__

    def __sub__(self, other):
        if isinstance(other, Pixel):
            x = self._x - other.x
            y = self._y - other.y
        elif isinstance(other, Iterable):
            return self - Pixel(*other)
        else:
            x = self._x - other
            y = self._y - other
        return Pixel(x,y)

    def __rsub__(self, other):
        if isinstance(other, Pixel):
            x = other.x - self._x
            y = other.y - self._y
        elif isinstance(other, Iterable):
            return Pixel(*other) - self
        else:
            x = other - self._x
            y = other - self._y
        return Pixel(x,y)

    def __mul__(self, other):
        if isinstance(other, Pixel):
            x = self._x * other.x
            y = self._y * other.y
        elif isinstance(other, Iterable):
            return self * Pixel(*other)
        else:
            x = self._x * other
            y = self._y * other
        return Pixel(x,y)

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, Pixel):
            x = self._x / other.x
            y = self._y / other.y
        elif isinstance(other, Iterable):
            return self / Pixel(*other)
        else:
            x = self._x / other
            y = self._y / other
        return Pixel(x,y)

    def __rtruediv__(self, other):
        if isinstance(other, Pixel):
            x = other.x / self._x
            y = other.y / self._y
        elif isinstance(other, Iterable):
            return Pixel(*other) / self
        else:
            x = other / self._x
            y = other / self._y
        return Pixel(x,y)

    def __floordiv__(self, other):
        if isinstance(other, Pixel):
            x = self._x // other.x
            y = self._y // other.y
        elif isinstance(other, Iterable):
            return self // Pixel(*other)
        else:
            x = self._x // other
            y = self._y // other
        return Pixel(x,y)

    def __rfloordiv__(self, other):
        if isinstance(other, Pixel):
            x = other.x // self._x
            y = other.y // self._y
        elif isinstance(other, Iterable):
            return Pixel(*other) // self
        else:
            x = other // self._x
            y = other // self._y
        return Pixel(x,y)

    __div__ = __truediv__
    __rdiv__ = __rtruediv__

    def as_tuple(self):
        return (self._x, self._y)

    def as_array(self):
        return np.array(self.as_tuple())

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
        return self._x

    @x.setter
    def x(self, x):
        self._x = x

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, y):
        self._y = y

#=============================================================================#
#===== Useful Functions ======================================================#
#=============================================================================#

def convert_pixels_to_units(value, pixel_size, magnification, units):
    """Converts a value or iterable from pixels to given units"""
    if units == 'pixels':
        return value
    elif units == 'microns':
        if isinstance(value, Iterable):
            return type(value)(np.array(value) * pixel_size / magnification)
        elif isinstance(value, Pixel):
            value.units = units
        return value * pixel_size / magnification
    else:
        raise RuntimeError('Incorrect string for units')

def convert_microns_to_units(value, pixel_size, magnification, units):
    """Converts a value or iterable from microns to given units"""
    if units == 'microns':
        return value
    elif units == 'pixels':
        if isinstance(value, Iterable):
            return type(value)(np.array(value) * magnification / pixel_size)
        elif isinstance(value, Pixel):
            value.units = units
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