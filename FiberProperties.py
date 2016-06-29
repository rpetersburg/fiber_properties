import numpy as np
from NumpyArrayHandler import NumpyArrayHandler

class FiberProperties(NumpyArrayHandler):

    def __init__(self, in_object=None, nf_object=None, ff_object=None):
        self.in_object = in_object
        self.nf_object = nf_object
        self.ff_object = ff_object
        self.nf_scrambling_gain = None
        self.ff_scrambling_gain = None

    def setInputObject(self, in_object):
        self.in_object = in_object

    def setNearFieldObject(self, nf_object):
        self.nf_object = nf_object

    def setFarFieldObject(self, ff_object):
        self.ff_object = ff_object

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
        in_center_y, in_center_x = in_object.getFiberCenterEdgeMethod()
        in_diameter = in_object.getFiberDiameter()

        out_centroid_y, out_centroid_x = out_object.getFiberCentroid()
        out_center_y, out_center_x = out_object.getFiberCenterEdgeMethod()
        out_diameter = out_object.getFiberDiameter()

        delta_D_in = np.sqrt((in_centroid_x - in_center_x)**2 + (in_centroid_y - in_center_y)**2)
        delta_D_out = np.sqrt((out_centroid_x - out_center_x)**2 + (out_centroid_y - out_center_y)**2)

        scramblingGain = (delta_D_in / in_diameter) / (delta_D_out / out_diameter)

        return scramblingGain

    def getModalNoise(self, radius_factor=0.95):
        nf_image = self.nf_object.getImageArray()

        radius = self.nf_object.getFiberRadius()
        center_y = self.nf_object.getFiberCenter()[0]
        center_x = self.nf_object.getFiberCenter()[1]

        intensity_list = []
        for i in xrange(self.nf_object.getImageWidth()):
            for j in xrange(self.nf_object.getImageHeight()):
                if (center_x-i)**2 + (center_y-j)**2 < (radius * radius_factor)**2:
                    intensity_list.append(nf_image[j,i])

        intensity_array = np.array(intensity_list)

        return intensity_array.var() / intensity_array.mean()
