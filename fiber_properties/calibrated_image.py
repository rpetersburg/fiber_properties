"""calibrated_image.py was written by Ryan Petersburg for use with fiber
characterization on the EXtreme PREcision Spectrograph
"""
import numpy as np
from .base_image import BaseImage
from .numpy_array_handler import filter_image, subframe_image

class CalibratedImage(BaseImage):
    """Fiber face image analysis class

    Class that contains calibration images and executes corrections based on
    those images

    Attributes
    ----------
    dark : str, array_like, or None
        The input used to set the dark image. See
        BaseImage.convert_image_to_array() for details
    ambient : str, array_like, or None
        The input used to set the ambient image. See
        BaseImage.convert_image_to_array() for details
    flat : str, array_like, or None
        The input used to set the flat image. See
        BaseImage.convert_image_to_array() for details
    kernel_size : int (odd)
        The kernel side length used when filtering the image. This value may
        need to be tweaked, especially with few co-added images, due to random
        noise. The filtered image is used for the centering algorithms, so for
        a "true test" use kernel_size=1, but be careful, because this may
        lead to needing a fairly high threshold for the noise.
    new_calibration : bool
        Whether or not self.calibration has been set with new images

    Args
    ----
    image_input : str, array_like, or None, optional
        See BaseImage class for details
    dark : str, array_like, or None, optional
        Image input to instantiate BaseImage for dark image
    ambient : str, array_like, or None, optional
        Image input to instantiate BaseImage for ambient image
    flat : str, array_like, or None, optional
        Image input to instantiate BaseImage for flat image
    kernel_size : int (odd), optional
        Set the kernel size for filtering
    **kwargs : keworded arguments
        Passed into the BaseImage superclass

    """
    def __init__(self, image_input, dark=None, ambient=None, flat=None,
                 kernel_size=9, **kwargs):
        self.dark = dark
        self.ambient = ambient
        self.flat = flat
        self.kernel_size = kernel_size
        self.new_calibration = True

        super(CalibratedImage, self).__init__(image_input, **kwargs)

    #=========================================================================#
    #==== Primary Image Getters ==============================================#
    #=========================================================================#

    def get_uncorrected_image(self):
        """Return the raw image without corrections or filtering.

        Returns
        -------
        uncorrected_image : 2D numpy array
            Raw image or average of images (depending on image_input)
        """
        return self.convert_image_to_array(self.image_input)

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
            return self.image_from_file(self.image_file)
        return self.execute_error_corrections(self.get_uncorrected_image())

    def get_uncorrected_filtered_image(self, kernel_size=None, **kwargs):
        """Return a median filtered image

        Args
        ----
        kernel_size : {None, int (odd)}, optional
            The side length of the kernel used to median filter the image. Uses
            self.kernel_size if None.

        Returns
        -------
        filtered_image : 2D numpy array
            The stored image median filtered with the given kernel_size
        """
        image = self.get_uncorrected_image()
        if image is None:
            return None
        if kernel_size is None:
            kernel_size = self.kernel_size
        return filter_image(image, kernel_size, **kwargs)

    def get_filtered_image(self, kernel_size=None, **kwargs):
        """Return an error corrected and median filtered image

        Returns
        -------
        filtered_image : 2D numpy array
            The stored image median filtered with the given kernel_size and
            error corrected using the given method
        """
        image = self.get_image()
        if image is None:
            return None
        if kernel_size is None:
            kernel_size = self.kernel_size
        return filter_image(image, kernel_size, **kwargs)

    #=========================================================================#
    #==== Calibration Image Getters ==========================================#
    #=========================================================================#

    def get_dark_image(self):
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
        return BaseImage(self.dark).get_image()

    def get_ambient_image(self):
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
        return CalibratedImage(self.ambient, dark=self.dark).get_image()

    def get_flat_image(self):
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
        return CalibratedImage(self.flat, dark=self.dark).get_image()

    def set_dark(self, dark):
        """Sets the dark calibration image."""
        self.dark = dark
        self.new_calibration = True

    def set_ambient(self, ambient):
        """Sets the ambient calibration image."""
        self.ambient = ambient
        self.new_calibration = True

    def set_flat(self, flat):
        """Sets the flat calibration images."""
        self.flat = flat
        self.new_calibration = True

    #=========================================================================#
    #==== Image Calibration Algorithm ========================================#
    #=========================================================================#

    def execute_error_corrections(self, image):
        """Applies corrective images to image

        Applies dark image to the flat field and ambient images. Then applies
        flat field and ambient image correction to the primary image

        Args
        ----
        image : 2D numpy array
            Image to be corrected

        Returns
        -------
        corrected_image : 2D numpy array
            Corrected image
        """
        if image is None:
            return None
        corrected_image = image

        dark_image = self.get_dark_image()
        if dark_image is not None and dark_image.shape != corrected_image.shape:
            dark_image = subframe_image(dark_image, self.subframe_x,
                                        self.subframe_y, self.width,
                                        self.height)
        corrected_image = self.remove_dark_image(corrected_image,
                                                 dark_image)

        ambient_image = self.get_ambient_image()
        if ambient_image is not None:
            if ambient_image.shape != corrected_image.shape:
                ambient_image = subframe_image(ambient_image, self.subframe_x,
                                               self.subframe_y, self.width,
                                               self.height)
            ambient_exp_time = BaseImage(self.ambient).exp_time
            if self.exp_time is not None and ambient_exp_time != self.exp_time:
                corrected_image = self.remove_dark_image(corrected_image,
                                                         ambient_image
                                                         * self.exp_time
                                                         / ambient_exp_time)
            else:
                corrected_image = self.remove_dark_image(corrected_image,
                                                         ambient_image)

        flat_image = self.get_flat_image()
        if flat_image is not None:
            if flat_image.shape != corrected_image.shape:
                flat_image = subframe_image(flat_image, self.subframe_x,
                                            self.subframe_y, self.width,
                                            self.height)
            corrected_image *= flat_image.mean() / flat_image

        self.new_calibration = False
        return corrected_image

    def remove_dark_image(self, image, dark_image=None):
        """Uses dark image to correct image

        Args
        ----
        image : 2D numpy array
            numpy array of the image
        dark_image : 2D numpy array
            dark image to be removed

        Returns
        -------
        output_array : 2D numpy array
            corrected image
        """
        if dark_image is None:
            dark_image = self.get_dark_image()
        if dark_image is None:
            dark_image = np.zeros_like(image)
        output_image = image - dark_image

        # Renormalize to the approximate smallest value (avoiding hot pixels)
        output_image -= filter_image(output_image, 5).min()
        # Prevent any dark/ambient image hot pixels from leaking through
        output_image *= (output_image > -1000.0).astype('uint8')

        return output_image

    #=========================================================================#
    #==== Attribute Setters ==================================================#
    #=========================================================================#

    def set_attributes_from_object(self, object_file):
        super(CalibratedImage, self).set_attributes_from_object(object_file)

        self.dark = self.change_path(self.dark)
        self.ambient = self.change_path(self.ambient)
        self.flat = self.change_path(self.flat)
