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
        dark
        flat
        ambient
    """
    def __init__(self, dark=None, ambient=None, flat=None):
        self.dark = dark
        self.ambient = ambient
        self.flat = flat

    def getDarkImage(self, full_output=False):
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
        return convertImageToArray(self.dark, full_output)

    def getAmbientImage(self, full_output=False):
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
        return convertImageToArray(self.ambient, full_output)

    def getFlatImage(self, full_output=False):
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
        return convertImageToArray(self.flat, full_output)

    def executeErrorCorrections(self, image, image_info=None,
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

        dark_image = self.getDarkImage()  
        if dark_image is None:
            dark_image = np.zeros_like(image)
        else:
            dark_image = subframeImage(dark_image, subframe_x, subframe_y,
                                       width, height)
        corrected_image = self.removeDarkImage(image, dark_image)

        ambient_image, ambient_info = self.getAmbientImage(True)   
        if ambient_image is not None:
            ambient_image = subframeImage(ambient_image, subframe_x,
                                          subframe_y, width, height)
            if 'exp_time' in ambient_info.__dict__:
                ambient_exp_time = ambient_info.exp_time 
            if exp_time is not None and ambient_exp_time is not None:
                corrected_image = self.removeDarkImage(corrected_image,
                                                       self.removeDarkImage(ambient_image,
                                                                            dark_image)
                                                       * exp_time / ambient_exp_time)
            else:
                corrected_image = self.removeDarkImage(corrected_image,
                                                       self.removeDarkImage(ambient_image,
                                                                            dark_image))

        flat_image = self.getFlatImage()
        if flat_image is not None:
            flat_image = subframeImage(flat_image, subframe_x,
                                       subframe_y, width, height)
            flat_image = self.removeDarkImage(flat_image, dark_image)
            corrected_image *= flat_image.mean() / flat_image

        return corrected_image

    def removeDarkImage(self, image_array, dark_image=None):
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
            dark_image = self.getDarkImage()
        output_array = image_array - dark_image

        # Prevent any pixels from becoming negative values
        output_array *= (output_array > 0.0).astype('float64')

        return output_array