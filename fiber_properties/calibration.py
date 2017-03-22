"""Calibration.py was written by Ryan Petersburg for use with fiber
characterization for the EXtreme PRecision Spectrograph
"""
import numpy as np
from fiber_properties.image_conversion import convert_image_to_array
from fiber_properties.numpy_array_handler import subframe_image
from fiber_properties.plotting import show_image_array

class Calibration(object):
    """Fiber face image analysis class

    Class that contains calibration images and executes corrections based on
        those images

    Attributes:
        dark
        flat
        ambient
    """
    def __init__(self, dark=None, ambient=None, flat=None):
        self.dark = dark
        self.ambient = ambient
        self.flat = flat

    def get_dark_image(self, full_output=False):
        """Returns the dark image.

        Args
        ----
        full_output : boolean, optional
            Passed to converImageToArray function

        Returns
        -------
        dark_image : 2D numpy array
            The dark image
        output_obj : ImageInfo, optional
            Object containing information about the image, if full_output=True
        """
        return convert_image_to_array(self.dark, full_output)

    def get_ambient_image(self, full_output=False):
        """Returns the ambient image.

        Args
        ----
        full_output : boolean, optional
            Passed to converImageToArray function

        Returns
        -------
        dark_image : 2D numpy array
            The dark image
        output_obj : ImageInfo, optional
            Object containing information about the image, if full_output=True
        """
        return convert_image_to_array(self.ambient, full_output)

    def get_flat_image(self, full_output=False):
        """Returns the flat image.

        Args
        ----
        full_output : boolean, optional
            Passed to converImageToArray function

        Returns
        -------
        dark_image : 2D numpy array
            The dark image
        output_obj : ImageInfo, optional
            Object containing information about the image, if full_output=True
        """
        return convert_image_to_array(self.flat, full_output)

    def execute_error_corrections(self, image, image_info=None,
                                  subframe_x=0, subframe_y=0, exp_time=None):
        """Applies corrective images to image

        Applies dark image to every instatiated image. Then applies flat field
        and ambient image correction to the primary image

        Args
        ----
        image : 2D numpy array
            Image to be corrected
        image_info : ImageInfo, optional
            ImageInfo object related to image

        Returns
        -------
        corrected_image : 2D numpy array
            Corrected image
        """
        height, width = image.shape
        if image_info is not None:
            subframe_x = image_info.subframe_x
            subframe_y = image_info.subframe_y
            exp_time = image_info.exp_time

        dark_image = self.get_dark_image()
        if dark_image is None:
            dark_image = np.zeros_like(image)
        elif dark_image.shape != image.shape:
            dark_image = subframe_image(dark_image, subframe_x, subframe_y,
                                        width, height)
        corrected_image = self.remove_dark_image(image, dark_image)

        ambient_image, ambient_info = self.get_ambient_image(True)
        if ambient_image is not None:
            if ambient_image.shape != image.shape:
                ambient_image = subframe_image(ambient_image, subframe_x,
                                               subframe_y, width, height)
            if 'exp_time' in ambient_info.__dict__ and exp_time is not None:
                ambient_exp_time = ambient_info.exp_time
                corrected_image = self.remove_dark_image(corrected_image,
                                                         self.remove_dark_image(ambient_image,
                                                                                dark_image)
                                                         * exp_time / ambient_exp_time)
            else:
                corrected_image = self.remove_dark_image(corrected_image,
                                                         self.remove_dark_image(ambient_image,
                                                                                dark_image))

        flat_image = self.get_flat_image()
        if flat_image is not None:
            flat_image = subframe_image(flat_image, subframe_x,
                                        subframe_y, width, height)
            flat_image = self.remove_dark_image(flat_image, dark_image)
            corrected_image *= flat_image.mean() / flat_image

        return corrected_image

    def remove_dark_image(self, image_array, dark_image=None):
        """Uses dark image to correct image

        Args
        ----
        image_array : 2D numpy array
            numpy array of the image

        Returns
        -------
        output_array : 2D numpy array
            corrected image
        """
        if dark_image is None:
            dark_image = self.get_dark_image()
        output_array = image_array - dark_image

        # Prevent any pixels from becoming negative values
        # output_array *= (output_array > 0.0).astype('float64')

        return output_array
