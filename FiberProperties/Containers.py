"""Containers.py was written by Ryan Petersburg for use with fiber
characterization on the EXtreme PREcision Spectrograph

The classes in this module are used as containers for information (similar to
dictionaries) in the ImageAnalysis class and ImageConcerversion functions.
These are used instead of dictionaries due to the simplicity of attribute
instantiation so that the information is ALWAYS either a value or NONE rather
than an empty slot in a dictionary.
"""
#=============================================================================#
#===== Metadata Containers ===================================================#
#=============================================================================#

class ImageInfo(object):
    """Container for an image's meta information"""
    def __init__(self):
        self.pixel_size = None
        self.camera = None
        self.magnification = None
        self.height = None
        self.width = None
        self.subframe_x = None
        self.subframe_y = None
        self.exp_time = None
        self.bit_depth = None
        self.date_time = None
        self.temp = None
        self.num_images = None
        self.folder = None
        self.test = None

class AnalysisInfo(object):
    """Container for meta information about ImageAnalysis."""
    def __init__(self, kernel_size, threshold):
        self.kernel_size = kernel_size
        self.threshold = threshold

class Edges(object):
    """Container for the fiber image edges."""
    def __init__(self):
        self.left = None
        self.right = None
        self.top = None
        self.bottom = None

class FiberInfo(object):
    """Container for information concerning the fiber grouped by method."""
    def __init__(self, info=None):
        if info == 'pixel':
            self.edge = Pixel()
            self.radius = Pixel()
            self.circle = Pixel()
            self.gaussian = Pixel()
            self.full = Pixel()
        elif info == 'value':
            self.edge = None
            self.radius = None
            self.circle = None
            self.gaussian = None
            self.full = None

class Pixel(object):
    """Container for the x and y position of a pixel."""
    def __init__(self, x=None, y=None):
        self.x = x
        self.y = y