"""fiber_image.py was written by Ryan Petersburg for use with fiber
characterization on the EXtreme PREcision Spectrograph
"""
import numpy as np
from .numpy_array_handler import (sum_array, isolate_circle, circle_array,
                                  polynomial_fit, gaussian_fit,
                                  rectangle_array, intensity_array,
                                  gaussian_array)
from .plotting import (plot_cross_sections, show_plots, plot_image)
from .containers import (FiberInfo, Edges, FRDInfo, ModalNoiseInfo,
                         convert_microns_to_units)
from .calibrated_image import CalibratedImage
from .modal_noise import modal_noise
from .fiber_centroid import fiber_centroid
from .fiber_center import fiber_center_and_diameter
from .fiber_edges import fiber_edges

#=============================================================================#
#===== FiberImage Class ======================================================#
#=============================================================================#

class FiberImage(CalibratedImage):
    """Fiber face image analysis class

    Class that conducts image analysis on a fiber face image after it has been
    corrected by the given dark, flat field, and ambient images. Also contains
    information about the CCD and camera that took the image. Public methods in
    this class allow calculation of the fiber face's centroid, center, and
    diameter using multiple different algorithms

    Attributes
    ----------
    _frd_info : FRDInfo
        Container for information pertaining to the calculated FRD
    _modal_noise_info : ModalNoiseInfo
        Container for information pertaining to the calculated modal noise

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

    gaussian_coeffs : tuple
        Coefficients of the gaussian fit function ordered as 
        (x0, y0, radius, amplitude, offset)

    Args
    ----
    image_input : str, array_like, or None
        See BaseImage class for details
    threshold : int, optional (default=256)
        The threshold value used with the centering method. This value may need
        to be tweaked for noisier images. Make sure to confirm decent results
        before setting this value too low
    input_fnum : float
        The focal ratio input into the FCS fiber
    output_fnum : float
        The focal ratio on the output side of the FCS
    **kwargs : keyworded arguments
        Passed into the CalibratedImage superclass
    """
    def __init__(self, image_input, threshold=1000,
                 input_fnum=2.4, output_fnum=2.4, **kwargs):
        # Private attribute initialization
        self.threshold = threshold

        self._frd_info = FRDInfo()
        self._frd_info.input_fnum = input_fnum
        self._frd_info.output_fnum = output_fnum

        self._modal_noise_info = ModalNoiseInfo()

        self.gaussian_coeffs = None

        self._edges = Edges()
        self._center = FiberInfo('pixel')
        self._centroid = FiberInfo('pixel')
        self._diameter = FiberInfo('value')

        super(FiberImage, self).__init__(image_input, **kwargs)

    #=========================================================================#
    #==== Fiber Data Getters =================================================#
    #=========================================================================#

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
        _center.method : .containers.Pixel
            The center of the fiber face in the context of the given method

        Returns
        -------
        _center.method : .containers.Pixel
            in the given units
        _diameter.method : float
            in the given units
        """
        center = self.get_fiber_center(method, units=units, **kwargs)
        kwargs['show_image'] = False
        diameter = self.get_fiber_diameter(method, units=units, **kwargs)
        return center, diameter

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
        _center.method : .containers.Pixel
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
        _center.method : .containers.Pixel
            The center of the fiber face in the context of the given method

        Returns
        -------
        _diameter.method : float
            in the given units
        """
        if method is None:
            if self.camera != 'ff':
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
        _center.method : .containers.Pixel
            The center of the fiber face in the context of the given method

        Returns
        -------
        _center.method : .containers.Pixel
            in the given units
        """
        if method is None:
            if self.camera != 'ff':
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

        if getattr(self._center, method).x is None or (method == 'circle' and
                                                       hasattr(kwargs, 'radius')):
            self.set_fiber_center(method, **kwargs)

        center = getattr(self._center, method)
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
        _centroid.method : .containers.Pixel
            The centroid of the fiber face in the context of the given method

        Returns
        -------
        _centroid.method : .containers.Pixel
            in the given units
        """
        if method is None:
            if self.camera != 'ff':
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

        centroid = getattr(self._centroid, method)
        return self.convert_pixels_to_units(centroid, units)

    #=========================================================================#
    #==== Image Fitting Getters ==============================================#
    #=========================================================================#

    def get_rectangle_fit(self):
        """Return the best rectangle fit for the image"""
        if self._center.circle.x is None:
            self.set_fiber_center(method='circle')
        rectangle_fit = rectangle_array(self.get_mesh_grid(),
                                        corners = [self._edges.left,
                                                   self._edges.bottom,
                                                   self._edges.right,
                                                   self._edges.top])
        return rectangle_fit

    def get_gaussian_fit(self, full_output=False, **kwargs):
        """Return the best gaussian fit for the image"""
        if self.gaussian_coeffs is None:
            self.set_gaussian_fit(**kwargs)

        mesh_grid = self.get_mesh_grid()
        coeffs = self.gaussian_coeffs
        gauss_fit = gaussian_array(mesh_grid, *coeffs).reshape(*mesh_grid[0].shape)

        if full_output:
            return gauss_fit, coeffs
        return gauss_fit

    def set_gaussian_fit(self, radius_factor=1.0, **kwargs):
        """Set the best gaussian fit coefficients for the image

        Args
        ----
        radius_factor : float (default=1.0)
            fraction of the fiber radius inside which the fit is made

        Sets
        ----
        self.gaussian_coeffs : tuple
            tuple of the gaussian coefficients ordered as 
            (x0, y0, radius, amplitude, offset)
        """
        image = self.get_image()
        center = self.get_fiber_center()
        radius = self.get_fiber_radius() * radius_factor
        if self.camera == 'in':
            initial_guess = (center.x, center.y, 100 / self.get_pixel_size(),
                             image.max(), image.min())
        else:
            initial_guess = (center.x, center.y, radius,
                             image.max(), image.min())

        gauss_fit, coeffs = gaussian_fit(image, initial_guess=initial_guess,
                                         full_output=True, center=center,
                                         radius=radius)

        self.gaussian_coeffs = coeffs

    def get_polynomial_fit(self, deg=6, radius_factor=0.95, **kwargs):
        """Return the best polynomial fit for the image

        Args
        ----
        deg : int (default=6)
            the degree of polynomial to fit
        radius_factor : float
            fraction of the fiber radius inside which the fit is shown

        Returns
        -------
        poly_fit : 2D numpy.ndarray
        """
        image = self.get_image()
        center = self.get_fiber_center(**kwargs)
        radius = self.get_fiber_radius(**kwargs) * radius_factor
        poly_fit = polynomial_fit(image, deg, center, radius)
        return poly_fit

    def get_tophat_fit(self, **kwargs):
        """Return the circle array that best covers the fiber face

        Args
        ----
        **kwargs :
            the keyworded arguments to pass to self.get_fiber_center()

        Returns
        -------
        circle_array : numpy.ndarray (2D)
            centered at best calculated center and with best calculated
            diameter
        """
        image = self.get_image()
        center = self.get_fiber_center(**kwargs)
        radius = self.get_fiber_radius(**kwargs)
        tophat_fit = circle_array(self.get_mesh_grid(), center.x, center.y,
                                  radius, res=1)
        tophat_fit *= np.median(intensity_array(image, center, radius))
        return tophat_fit

    #=========================================================================#
    #==== Focal Ratio Degradation Methods ====================================#
    #=========================================================================#

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
            Container for FRD information. See containers.FRDInfo
        """
        if new or not self._frd_info.energy_loss:
            self.set_frd_info(**kwargs)
        return self._frd_info

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
            the approximate focal ratio inside which fnum_diameter of the total
            encircled energy is included
        _frd_info.energy_loss : float
            the loss of energy when the output focal ratio equals the input
            focal ratio given as a percent
        _frd_info.encircled_energy_fnum : list(float)
            list of the focal ratios used to calculate encircled energy
        _frd_info.encircled_energy : list(float)
            list of the encircled energy at each given focal ratio
        """
        center = self.get_fiber_centroid(method='full')

        fnums = list(np.arange(f_lim[0], f_lim[1] + res, res))
        energy_loss = None
        output_fnum = None
        encircled_energy = []
        image = self.get_image()
        for fnum in fnums:
            radius = self.convert_fnum_to_radius(fnum, units='pixels')
            isolated_circle = isolate_circle(image,
                                             center,
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
    #==== Modal Noise Methods ================================================#
    #=========================================================================#

    def get_modal_noise(self, method='fft', new=False, **kwargs):
        if (not hasattr(self._modal_noise_info, method) 
            or getattr(self._modal_noise_info, method) is None) or new:
            self.set_modal_noise(method, **kwargs)
        return getattr(self._modal_noise_info, method)

    def set_modal_noise(self, method=None, **kwargs):
        if method is None:
            if self.camera == 'nf':
                method1 = 'tophat'
            elif self.camera == 'ff':
                method1 = 'gaussian'
            methods = [method1, 'polynomial', 'contrast', 'filter', 'gradient', 'fft']
        else:
            methods = [method]

        for method in methods:
            setattr(self._modal_noise_info, method,
                    modal_noise(self, method, **kwargs))            

    #=========================================================================#
    #==== Image Centroiding ==================================================#
    #=========================================================================#

    def set_fiber_centroid(self, method='full', **kwargs):
        """Find the centroid of the fiber face image

        Args
        ----
        method : {'full', 'edge', 'radius', 'gaussian', 'circle'}, optional
            If 'full', takes the centroid of the entire image. Otherwise, uses
            the specified method to isolate only the fiber face in the image
        radius_factor : number, optional
            The factor by which the radius is multiplied when isolating the
            fiber face in the image
        show_image : bool, optional
            Shows centroid dot on fiber image
        fiber_shape : {'circle', 'rectangle'}, optional
            The shape of the fiber core cross-section. Used to decide which
            points to use when calculating the centroid.

        Sets
        ----
        _centroid.method : .containers.Pixel
            The centroid of the image in the context of the given method
        """
        setattr(self._centroid, method, fiber_centroid(self, method, **kwargs))

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
        _centroid.method : .containers.Pixel
            The centroid of the image in the context of the given method
        _center.method : .containers.Pixel
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
        _center.method : .containers.Pixel
            The center of the fiber face in the context of the given method
        """
        if method == 'circle':
            raise RuntimeError('Fiber diameter cannot be set by circle method')
        self.set_fiber_center(method, **kwargs)

    def set_fiber_center(self, method, **kwargs):
        """Find fiber center using given method

        Args
        ----
        method : {'edge', 'radius', 'gaussian', 'circle'}
            Uses the respective method to find the fiber center
        **kwargs :
            The keyworded arguments to pass to the centering method

        Sets
        ----
        _center.method : .containers.Pixel
            The center of the fiber face in the context of the given method
        _diameter.method : float
            The diameter of the fiber face in the context of the given method
        """
        center, diameter = fiber_center_and_diameter(self, method, **kwargs)
        setattr(self._center, method, center)
        setattr(self._diameter, method, diameter)

    def set_fiber_edges(self, **kwargs):
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
        """
        self._edges = fiber_edges(self, **kwargs)
        # image = self.get_filtered_image() # To prvent hot pixels

        # left = -1
        # right = -1
        # for index in xrange(self.width):
        #     if left < 0:
        #         if image[:, index].max() > self.threshold:
        #             left = index
        #     else:
        #         if image[:, index].max() > self.threshold:
        #             right = index

        # top = -1
        # bottom = -1
        # for index in xrange(self.height):
        #     if top < 0:
        #         if image[index, :].max() > self.threshold:
        #             top = index
        #     else:
        #         if image[index, :].max() > self.threshold:
        #             bottom = index

        # left = np.array([left, image[:, left].argmax()])
        # right = np.array([right, image[:, right].argmax()])
        # top = np.array([image[top, :].argmax(), top])
        # bottom = np.array([image[bottom, :].argmax(), bottom])
        # # diameter = (np.sqrt(((right - left)**2).sum())
        # #           + np.sqrt(((bottom - top)**2).sum())) / 2.0

        # self._edges.left.set_pixel(*left)
        # self._edges.right.set_pixel(*right)
        # self._edges.top.set_pixel(*top)
        # self._edges.bottom.set_pixel(*bottom)
        # # self._diameter.edge = diameter

    #=========================================================================#
    #==== Useful Methods =====================================================#
    #=========================================================================#

    def convert_fnum_to_radius(self, fnum, units):
        """Returns the value in the proper units"""
        return convert_fnum_to_radius(fnum,
                                      self.pixel_size,
                                      self.magnification,
                                      units)

    def plot_cross_sections(self, image=None, row=None, column=None):
        """Plots cross sections across the center of the fiber"""
        if image is None:
            image = self.get_image()
        if row is None:
            row = self.get_fiber_center().y
        if column is None:
            column = self.get_fiber_center().x
        plot_cross_sections(image, row, column)

#=============================================================================#
#===== Useful Functions ======================================================#
#=============================================================================#

def convert_fnum_to_radius(fnum, pixel_size, magnification, units='pixels'):
    """Converts a focal ratio to an image radius in given units."""
    fcs_focal_length = 4.0 # inches
    diameter = fcs_focal_length / fnum # inches
    radius = 25400 * diameter / 2.0 # microns
    return convert_microns_to_units(radius, pixel_size, magnification, units)

if __name__ == "__main__":
    folder = 'C:/Libraries/Box Sync/ExoLab/Fiber_Characterization/Image Analysis/data/scrambling/2016-08-05 Prototype Core Extension 1/'

    images = [folder + 'Shift_00/in_' + str(i).zfill(3) + '.fit' for i in xrange(10)]
    dark = [folder + 'Dark/in_' + str(i).zfill(3) + '.fit' for i in xrange(10)]
    ambient = [folder + 'Ambient/in_' + str(i).zfill(3) + '.fit' for i in xrange(10)]

    im_obj = FiberImage(images, dark=dark, ambient=ambient)

    tol = 0.25
    test_range = 5
    factor = 1.0

    im_obj.show_image()
    print
    print 'Centroid:'
    print im_obj.get_fiber_centroid(method='full', radius_factor=factor)
    print
    print 'Edge:'
    print 'center:', im_obj.get_fiber_center(method='edge')
    print 'centroid:', im_obj.get_fiber_centroid(method='edge', radius_factor=factor)
    print 'diameter:', im_obj.get_fiber_diameter(method='edge', units='microns'), 'microns'
    print
    print 'Radius:'
    print 'center:', im_obj.get_fiber_center(method='radius', tol=tol, test_range=test_range)
    print 'centroid:', im_obj.get_fiber_centroid(method='radius', radius_factor=factor)
    print 'diameter:', im_obj.get_fiber_diameter(method='radius', units='microns'), 'microns'
    print
    print 'Gaussian:'
    print 'center:', im_obj.get_fiber_center(method='gaussian')
    print 'centroid:', im_obj.get_fiber_centroid(method='gaussian', radius_factor=factor)
    print 'diameter:', im_obj.get_fiber_diameter(method='gaussian', units='microns'), 'microns'

    im_obj.save()

    new_im_obj = FiberImage(im_obj.object_file)

    new_im_obj.show_image()
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
