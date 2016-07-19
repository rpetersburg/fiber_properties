import numpy as np
from NumpyArrayHandler import NumpyArrayHandler

class ScramblingGain(NumpyArrayHandler):

    def __init__(self, in_object=None, nf_object=None, ff_object=None):
        self.in_object = in_object
        self.nf_object = nf_object
        self.ff_object = ff_object
        self.nf_scrambling_gain = None
        self.ff_scrambling_gain = None

#=============================================================================#
#==== Private Variable Setters ===============================================#
#=============================================================================#

    def setInputObject(self, in_object):
        self.in_object = in_object

    def setNearFieldObject(self, nf_object):
        self.nf_object = nf_object

    def setFarFieldObject(self, ff_object):
        self.ff_object = ff_object

#=============================================================================#
#==== Scrambling Gain Methods ================================================#
#=============================================================================#

    def getNearFieldScramblingGain(self):
        if self.nf_scrambling_gain is None:
            self.nf_scrambling_gain = getScramblingGain(self.in_object, self.nf_object)
        return self.nf_scrambling_gain

    def getFarFieldScramblingGain(self):
        if self.ff_scrambling_gain is None:
            self.ff_scrambling_gain = getScramblingGain(self.in_object, self.ff_object)
        return self.ff_scrambling_gain

    def getScramblingGain(self, in_object, out_object):
        in_centroid_y, in_centroid_x = in_object.getFiberCentroid()
        in_y0, in_x0 = in_object.getFiberCenterEdgeMethod()
        in_diameter = in_object.getFiberDiameter()

        out_centroid_y, out_centroid_x = out_object.getFiberCentroid()
        out_y0, out_x0 = out_object.getFiberCenterEdgeMethod()
        out_diameter = out_object.getFiberDiameter()

        delta_D_in = np.sqrt((in_centroid_x - in_x0)**2 + (in_centroid_y - in_y0)**2)
        delta_D_out = np.sqrt((out_centroid_x - out_x0)**2 + (out_centroid_y - out_y0)**2)

        scramblingGain = (delta_D_in / in_diameter) / (delta_D_out / out_diameter)

        return scramblingGain

#=============================================================================#
#==== Show Methods ===========================================================#
#=============================================================================#

    def showInputImageArray(self):
        self.showImageArray(self.in_object.getImageArray())

    def showNearFieldImageArray(self):
        self.showImageArray(self.nf_object.getImageArray())

    def showFarFieldImageArray(self):
        self.showImageArray(self.ff_object.getImageArray())

    def showImageArray(self, image_array=None):
        if image_array is not None:
            super(FiberProperties, self).showImageArray(image_array)
        else:
            if self.in_object is not None:
                self.showInputImageArray()
            if self.nf_object is not None:
                self.showNearFieldImageArray()
            if self.ff_object is not None:
                self.showFarFieldImageArray()
