"""Calibration.py was written by Ryan Petersburg for use with fiber
characterization for the EXtreme PRecision Spectrograph
"""
import numpy as np
from ImageConversion import convertImageToArray
from NumpyArrayHandler import subframeImage

class Calibration(object):
    """Fiber face image analysis class

    Class that contains calibration images and executes corrections based on
        those images

    Attributes:
        dark_image
        flat_image
        ambient_image
    """
    def __init__(self, dark=None, flat=None, ambient=None):
        self.setDarkImage(dark)
        self.setFlatImage(flat)
        self.setAmbientImage(ambient)

    def setDarkImage(self, image_input):
        """Sets the corrective dark image

        Args:
            image_input: see convertImageToArray for options

        Sets:
            self.dark_image
        """
        self.dark_image = convertImageToArray(image_input)

    def setFlatImage(self, image_input):
        """Sets the corrective flat field image

        Args:
            image_input: see convertImageToArray for options

        Sets:
            self.flat_image
        """
        self.flat_image = convertImageToArray(image_input)

    def setAmbientImage(self, image_input):
        """Sets the corrective ambient image

        Args:
            *image_input: see convertImageToArray for options

        Sets:
            self.ambient_image
        """
        self.ambient_image, image_info = convertImageToArray(image_input, full_output=True)
        if 'exp_time' in image_info.__dict__:
            self._ambient_exp_time = image_info.exp_time

    def executeErrorCorrections(self, image, subframe_x=0, subframe_y=0, exp_time=None):
        """Applies corrective images to image

        Applies dark image to every instatiated image. Then applies flat field
        and ambient image correction to the primary image
        """
        height, width = image.shape

        if self.dark_image is None:
            temp_dark_image = np.zeros_like(image)
        else:
            temp_dark_image = subframeImage(self.dark_image, subframe_x,
                                            subframe_y, width, height)
        corrected_image = self.removeDarkImage(image, temp_dark_image)

        if self.ambient_image is not None:
            temp_ambient_image = subframeImage(self.ambient_image, subframe_x,
                                               subframe_y, width, height)
            if exp_time is None:
                corrected_image = self.removeDarkImage(corrected_image,
                                                       self.removeDarkImage(temp_ambient_image,
                                                                            temp_dark_image))
            else:
                corrected_image = self.removeDarkImage(corrected_image,
                                                       self.removeDarkImage(temp_ambient_image,
                                                                            temp_dark_image)
                                                       * exp_time / self._ambient_exp_time)

        if self.flat_image is not None:
            temp_flat_image = subframeImage(self.flat_image, subframe_x,
                                            subframe_y, width, height)
            corrected_flat_image = self.removeDarkImage(self.flat_image, temp_dark_image)
            corrected_image *= corrected_flat_image.mean() / corrected_flat_image

        return corrected_image

    def removeDarkImage(self, image_array, dark_image=None):
        """Uses dark image to correct image

        Args:
            image_array: numpy array of the image

        Returns:
            image_array: corrected image
        """
        if dark_image is None:
            dark_image = self.dark_image
        output_array = image_array - dark_image

        # Prevent any pixels from becoming negative values
        output_array *= (output_array > 0.0).astype('float64')

        return output_array