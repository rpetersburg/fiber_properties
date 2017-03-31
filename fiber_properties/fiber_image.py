"""fiber_image.py was written by Ryan Petersburg for use with fiber
characterization on the EXtreme PREcision Spectrograph
"""
import numpy as np
from .numpy_array_handler import (sum_array, crop_image, remove_circle,
                                  isolate_circle, gaussian_array, circle_array,
                                  polynomial_fit, gaussian_fit, rectangle_array)
from .plotting import (plot_cross_sections, plot_overlaid_cross_sections,
                       plot_dot, show_plots, plot_image_array)
from .containers import FiberInfo, Edges, FRDInfo
from .calibrated_image import CalibratedImage
from .base_image import convert_microns_to_units

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

    _gaussian_amp : float
        Amplitude of the gaussian fit function
    _gaussian_offset : float
        Offset of the gaussian fit function

    _phi : float
        Golden ratio for the optimization tests

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
    def __init__(self, image_input, threshold=256,
                 input_fnum=2.4, output_fnum=2.4, **kwargs):
        # Private attribute initialization
        super(FiberImage, self).__init__(image_input, **kwargs)

        self.threshold = threshold

        self._edges = Edges()
        self._center = FiberInfo('pixel')
        self._centroid = FiberInfo('pixel')
        self._diameter = FiberInfo('value')
        self._array_sum = FiberInfo('value')

        self._frd_info = FRDInfo()
        self._frd_info.input_fnum = input_fnum
        self._frd_info.output_fnum = output_fnum

        self._phi = (5 ** 0.5 - 1) / 2

        self._gaussian_amp = 0.0
        self._gaussian_offset = 0.0

        self._rectangle_width = 0.0
        self._rectangle_height = 0.0
        self._rectangle_angle = 0.0

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
        _center.method : {'x': float, 'y': float}
            The center of the fiber face in the context of the given method

        Returns
        -------
        _center.method.y : float
            in the given units
        _center.method.x : float
            in the given units
        _diameter.method : float
            in the given units
        """
        center_y, center_x = self.get_fiber_center(method, units=units, **kwargs)
        kwargs['show_image'] = False
        diameter = self.get_fiber_diameter(method, units=units, **kwargs)
        return center_y, center_x, diameter

    def get_fiber_radius(self, method=None, units='pixels', **kwargs):
        """Return the fiber radius

        Finds the radius of the fiber using the given method or, if no method
        is given, the most precise method already completed

        Args
        ----
        method : None or str {'radius', 'gaussian', 'circle', 'edge'}, optional
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
        _center.method : {'x': float, 'y': float}
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
        _center.method : {'x': float, 'y': float}
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
        _center.method : {'x': float, 'y': float}
            The center of the fiber face in the context of the given method

        Returns
        -------
        _center.method.y : float
            in the given units
        _center.method.x : float
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

        if getattr(self._center, method).x is None or method == 'circle':
            self.set_fiber_center(method, **kwargs)

        center = (getattr(self._center, method).y,
                  getattr(self._center, method).x)

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
        _centroid.method : {'x': float, 'y': float}
            The centroid of the fiber face in the context of the given method

        Returns
        -------
        _centroid.method.y : float
            in the given units
        _centroid.method.x : float
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
        centroid = (getattr(self._centroid, method).y,
                    getattr(self._centroid, method).x)

        return self.convert_pixels_to_units(centroid, units)

    #=========================================================================#
    #==== Image Fitting Getters ==============================================#
    #=========================================================================#

    def get_rectangle_fit(self):
        """Return the best rectangle fit for the image"""
        if self._center.rectangle.x is None:
            self.set_fiber_center(method='rectangle')
        rectangle_fit = rectangle_array(self.get_mesh_grid(),
                                        self._center.rectangle.x,
                                        self._center.rectangle.y,
                                        self._rectangle_width,
                                        self._rectangle_height,
                                        self._rectangle_angle
                                       ).reshape(self.get_height(),
                                                 self.get_width())
        return rectangle_fit

    def get_gaussian_fit(self):
        """Return the best gaussian fit for the image

        Returns
        -------
        _fit.gaussian : 2D numpy.ndarray
        """
        if self._center.gaussian.x is None:
            self.set_fiber_center(method='gaussian')
        gaussian_fit = gaussian_array(self.get_mesh_grid(),
                                      self._center.gaussian.x,
                                      self._center.gaussian.y,
                                      self._diameter.gaussian / 2.0,
                                      self._gaussian_amp,
                                      self._gaussian_offset
                                     ).reshape(self.get_height(),
                                               self.get_width())

        if self.camera == 'in':
            y0, x0 = self.get_fiber_center()
            radius = self.get_fiber_radius()
            filtered_image = self.get_filtered_image()
            fiber_face = circle_array(self.get_mesh_grid(), x0, y0, radius)
            fiber_face *= np.median(crop_image(filtered_image, x0, y0, radius)[0])
            gaussian_fit += fiber_face

        return gaussian_fit

    def get_polynomial_fit(self, deg=6, x0=None, y0=None):
        """Return the best polynomial fit for the image

        Args
        ----
        deg : int (default=6)
            The degree of polynomial to fit
        x0 : number
            The center column to use for the radial polynomial. Uses best
            calculated center if None.
        y0 : number
            The center row to use for the radial polynomial. Uses best
            calculated center if None.

        Returns
        -------
        polynomial_fit : 2D numpy.ndarray
        """
        if y0 is None or x0 is None:
            y0, x0 = self.get_fiber_center()
        return polynomial_fit(self.get_image(), deg, x0, y0)

    def get_tophat_fit(self):
        """Return the circle array that best covers the fiber face

        Returns
        -------
        circle_array : numpy.ndarray (2D)
            Circle array centered at best calculated center and with best
            calculated diameter
        """
        y0, x0 = self.get_fiber_center()
        radius = self.get_fiber_radius()
        return self.get_image().max() * circle_array(self.get_mesh_grid(), x0, y0, radius, res=1)

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
            Container for FRD information. See Containers.FRDInfo
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
            the approximate focal ratio inside which 95% of the total encircled energy
            is included
        _frd_info.energy_loss : float
            the loss of energy when the output focal ratio equals the input focal ratio
            given as a percent
        _frd_info.encircled_energy_fnum : list(float)
            list of the focal ratios used to calculate encircled energy
        _frd_info.encircled_energy : list(float)
            list of the encircled energy at each given focal ratio
        """
        center_y, center_x = self.get_fiber_centroid(method='full')

        fnums = list(np.arange(f_lim[0], f_lim[1] + res, res))
        energy_loss = None
        output_fnum = None
        encircled_energy = []
        image = self.get_image()
        for fnum in fnums:
            radius = self.convert_fnum_to_radius(fnum, units='pixels')
            isolated_circle = isolate_circle(image,
                                             center_x,
                                             center_y,
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
    #==== Image Centroiding ==================================================#
    #=========================================================================#

    def set_fiber_centroid(self, method='full', radius_factor=1.0,
                           show_image=False, **kwargs):
        """Find the centroid of the fiber face image

        Args
        ----
        method : {'full', 'edge', 'radius', 'gaussian', 'circle'}, optional
            If 'full', takes the centroid of the entire image. Otherwise, uses
            the specified method to isolate only the fiber face in the image
        radius_factor : number
            The factor by which the radius is multiplied when isolating the
            fiber face in the image

        Sets
        ----
        _centroid.method : Pixel
            The centroid of the image in the context of the given method
        """
        image = self.get_filtered_image()
        if method == 'full':
            image_array_iso = image
        else:
            y0, x0 = self.get_fiber_center(method=method,
                                           show_image=False,
                                           **kwargs)
            radius = self.get_fiber_radius(method=method,
                                           show_image=False,
                                           **kwargs)
            image_array_iso = isolate_circle(image, x0, y0,
                                             radius*radius_factor, res=1)

        x_array, y_array = self.get_mesh_grid()
        getattr(self._centroid, method).x = ((image_array_iso * x_array).sum()
                                             / image_array_iso.sum())
        getattr(self._centroid, method).y = ((image_array_iso * y_array).sum()
                                             / image_array_iso.sum())

        if show_image:
            if method == 'gaussian':
                plot_overlaid_cross_sections(image, self.get_gaussian_fit(),
                                             self._center.gaussian.y,
                                             self._center.gaussian.x)
            plot_dot(image,
                     getattr(self._centroid, method).y,
                     getattr(self._centroid, method).x)
            show_plots()

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
        _centroid.method : Pixel
            The centroid of the image in the context of the given method
        _center.method : Pixel
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
        _center.method : Pixel
            The center of the fiber face in the context of the given method

        Raises
        ------
        RuntimeError
            cannot accept the 'circle' method when setting the diameter since
            it requires a known radius to run
        """
        if method == 'circle':
            raise RuntimeError('Fiber diameter cannot be set by circle method')
        self.set_fiber_center(method, **kwargs)

    def set_fiber_center(self, method, show_image=False, **kwargs):
        """Find fiber center using given method

        Args
        ----
        method : {'edge', 'radius', 'gaussian', 'circle'}
            Uses the respective method to find the fiber center
        show_image : boolean, optional (default=False)
            Whether or not to show relevant fitting image
        **kwargs :
            The keyworded arguments to pass to the centering method

        Raises
        ------
        RuntimeError
            needs a valid method string to run the proper algorithm
        """
        # Reset the fits due to new fiber parameters
        if method == 'radius':
            self.set_fiber_center_radius_method(**kwargs)
        elif method == 'edge':
            self.set_fiber_center_edge_method()
        elif method == 'circle':
            self.set_fiber_center_circle_method(**kwargs)
        elif method == 'gaussian':
            self.set_fiber_center_gaussian_method()
        elif method == 'rectangle':
            self.set_fiber_center_rectangle_method(**kwargs)
        else:
            raise RuntimeError('Incorrect string for fiber centering method')

        if show_image:
            x0 = getattr(self._center, method).x
            y0 = getattr(self._center, method).y
            r = getattr(self._diameter, method) / 2.0
            image = self.get_filtered_image()

            if method == 'gaussian':
                plot_overlaid_cross_sections(image, self.get_gaussian_fit(),
                                             y0, x0)
                plot_dot(image, y0, x0)
                show_plots()
            else:
                plot_image_array(remove_circle(image, x0, y0, r, res=1))
                plot_overlaid_cross_sections(image, image.max() / 2.0
                                             * circle_array(self.get_mesh_grid(),
                                                            x0, y0, r, res=1),
                                             y0, x0)
                show_plots()

    def set_fiber_center_rectangle_method(self, radius=None, **kwargs):
        """Set fiber center using a rectangle mask

        Uses Scipy.optimize.curve_fit method to fit fiber image to
        rectangle_array().

        Sets
        ----
        _diameter.rectangle : float
            sqrt(height^2 + width^2) of the rectangle)
        _center.rectangle : Pixel()
            Center of the fiber in the rectangle method context
        _fit.rectangle : 2D numpy.ndarray
            Best rectangle fit for the fiber image
        """
        if self._edges.left is None:
            self.set_fiber_edges()

        image = self.get_filtered_image()

        left = np.array([image[:, self._edges.left].argmax(),
                         self._edges.left])
        right = np.array([image[:, self._edges.right].argmax(),
                          self._edges.right])
        top = np.array([self._edges.top,
                        image[self._edges.top, :].argmax()])
        bottom = np.array([self._edges.bottom,
                           image[self._edges.bottom, :].argmax()])

        radius = (np.sqrt(((right - left)**2).sum())
                  + np.sqrt(((bottom - top)**2).sum())) / 4.0

        self.set_fiber_center_circle_method(radius, **kwargs)

        self._center.rectangle.x = self._center.circle.x
        self._center.rectangle.y = self._center.circle.y
        self._diameter.rectangle = radius * 2.0

    def set_fiber_center_gaussian_method(self):
        """Set fiber center using a Gaussian Fit

        Uses Scipy.optimize.curve_fit method to fit fiber image to
        gaussian_array(). The radius found extends to 2-sigma of the gaussian
        therefore encompassing ~95% of the imaged light. Use previous methods
        of center-finding to approximate the location of the center

        Sets
        ----
        _diameter.gaussian : float
            Diameter of the fiber in the gaussian method context
        _center.gaussian : {'x': float, 'y': float}
            Center of the fiber in the gaussian method context
        _fit.gaussian : 2D numpy.ndarray
            Best gaussian fit for the fiber image
        """
        fiber_y0, fiber_x0 = self.get_fiber_center(method='edge')
        fiber_radius = self.get_fiber_radius(method='edge')
        filtered_image = self.get_filtered_image()

        if self.camera == 'in':
            radius = -1
            factor = 0.9
            fiber_face = circle_array(self.get_mesh_grid(), fiber_x0, fiber_y0,
                                      fiber_radius, res=1)
            fiber_face *= np.median(crop_image(filtered_image, fiber_x0,
                                               fiber_y0, fiber_radius)[0])
            while radius < 1 or radius > fiber_radius:
                factor += 0.1
                cropped_image, new_x0, new_y0 = crop_image(filtered_image,
                                                           fiber_x0, fiber_y0,
                                                           fiber_radius*factor)

                initial_guess = (new_x0, new_y0, 100 / self.get_pixel_size(),
                                 cropped_image.max(), cropped_image.min())
                try:
                    fit, opt_parameters = gaussian_fit(cropped_image,
                                                       initial_guess=initial_guess,
                                                       full_output=True)

                    radius = abs(opt_parameters[2])
                except RuntimeError:
                    radius = -1

            x0 = opt_parameters[0] + int(fiber_x0-fiber_radius*factor)
            y0 = opt_parameters[1] + int(fiber_y0-fiber_radius*factor)
            amp = opt_parameters[3]
            offset = opt_parameters[4]

        else:
            initial_guess = (fiber_x0, fiber_y0, fiber_radius,
                             filtered_image.max(), filtered_image.min())

            _, opt_parameters = gaussian_fit(filtered_image,
                                             initial_guess=initial_guess,
                                             full_output=True)
            x0 = opt_parameters[0]
            y0 = opt_parameters[1]
            radius = abs(opt_parameters[2])
            amp = opt_parameters[3]
            offset = opt_parameters[4]

        self._center.gaussian.x = x0
        self._center.gaussian.y = y0
        self._diameter.gaussian = radius * 2.0
        self._gaussian_amp = amp
        self._gaussian_offset = offset

    def set_fiber_center_radius_method(self, radius_tol=.03, radius_range=None, **kwargs):
        """Set fiber center using dark circle with varying radius

        Uses a golden mean optimization method to find the optimal radius of the
        dark circle that covers the fiber image used in
        get_fiber_centerCircleMethod(). The optimization is for a parameter
        array_sum which is weighted by the area of the circle, meaning that a
        smaller circle is preferred over one that simply covers the entire image

        Args
        ----
        radius_tol : number (default=1)
            Minimum possible range of radius values before ending iteration
        radius_range: int (in pixels)
            Range of tested radii, i.e. max(radius) - min(radius). If None,
            uses full possible range

        Sets
        ----
        _diameter.radius : float
            Diameter of the fiber in the radius method context
        _center.radius : {'x': float, 'y': float}
            Center of the fiber in the radius method context
        _diameter.circle : float
            Also uses the circle method, therefore changes this value
        _center.circle : float
            Also uses the circle method, therefore chnages this value
        """
        image = self.get_filtered_image()

        # Initialize range of tested radii
        r = np.zeros(4).astype(float)

        if radius_range is not None:
            approx_radius = self.get_fiber_radius()
            radius_range /= 2.0

            r[0] = approx_radius - radius_range
            if r[0] < 0.0:
                r[0] = 0.0
            r[3] = approx_radius + radius_range
        else:
            r[0] = 0
            r[3] = min(self.height, self.width) / 2.0

        r[1] = r[0] + (1 - self._phi) * (r[3] - r[0])
        r[2] = r[0] + self._phi * (r[3] - r[0])

        array_sum = np.zeros(2).astype(float)
        for i in xrange(2):
            self.set_fiber_center(method='circle', radius=r[i+1],
                                  image=image, **kwargs)
            array_sum[i] = (self._array_sum.circle
                            + self.threshold
                            * np.pi * r[i+1]**2)

        min_index = np.argmin(array_sum) # Integer 0 or 1 for min of r[1], r[2]

        while abs(r[3]-r[0]) > radius_tol:
            if min_index == 0:
                r[3] = r[2]
                r[2] = r[1]
                r[1] = r[0] + (1 - self._phi) * (r[3] - r[0])
            else:
                r[0] = r[1]
                r[1] = r[2]
                r[2] = r[0] + self._phi * (r[3] - r[0])

            array_sum[1 - min_index] = array_sum[min_index]

            self.set_fiber_center(method='circle', radius=r[min_index+1],
                                  image=image, **kwargs)
            array_sum[min_index] = (self._array_sum.circle
                                    + self.threshold
                                    * np.pi * r[min_index+1]**2)

            min_index = np.argmin(array_sum) # Integer 0 or 1 for min of r[1], r[2]

        self._diameter.radius = r[min_index+1] * 2
        self._center.radius.y = self._center.circle.y
        self._center.radius.x = self._center.circle.x
        self._array_sum.radius = np.amin(array_sum)

    def set_fiber_center_circle_method(self, radius=None, center_tol=.03,
                                       center_range=None, image=None):
        """Finds fiber center using a dark circle of set radius

        Uses golden mean method to find the optimal center for a circle
        covering the fiber image. The optimization is for a parameter array_sum
        that simply sums over the entire fiber image array

        Args
        ----
        radius : float
            Radius to use when creating circle
        center_tol : number (default=1)
            Minimum possible range of center values before ending iteration
        center_range: int (in pixels)
            Range of tested centers, i.e. max(x0) - min(x0). If None,
            uses full possible range
        image : 2d numpy.ndarray, optional
            The image being analyzed. This is only useful for the radius_method.
            Probably not for use outside the class.

        Sets
        ----
        _diameter.circle : float
            Diameter of the fiber in the circle method context
        _center.circle : {'x': float, 'y': float}
            Center of the fiber in the circle method context
        _diameter.edge : float
            If center_range is not None, approximates the circle's center using
            the edge method
        _center.edge : float
            If center_range is not None, approximates the circle's center using
            the edge method
        """
        res = int(1.0/center_tol)
        if image is None:
            image = self.get_filtered_image()
        if radius is None:
            if self._center.circle.x is not None:
                return
            radius = self.get_fiber_radius(method='edge')
        print radius

        # Create four "corners" to test center of the removed circle
        x = np.zeros(4).astype(float)
        y = np.zeros(4).astype(float)

        if center_range is not None:
            approx_center = self.get_fiber_center(method='edge', show_image=False)
            center_range = center_range / 2.0

            x[0] = approx_center[1] - center_range
            if x[0] < radius:
                x[0] = radius
            x[3] = approx_center[1] + center_range
            if x[3] > self.width - radius:
                x[3] = self.width - radius

            y[0] = approx_center[0] - center_range
            if y[0] < radius:
                y[0] = radius
            y[3] = approx_center[0] + center_range
            if y[3] > self.height - radius:
                y[3] = self.height - radius

        else:
            x[0] = radius
            x[3] = self.width - radius

            y[0] = radius
            y[3] = self.height - radius

        x[1] = x[0] + (1 - self._phi) * (x[3] - x[0])
        x[2] = x[0] + self._phi * (x[3] - x[0])

        y[1] = y[0] + (1 - self._phi) * (y[3] - y[0])
        y[2] = y[0] + self._phi * (y[3] - y[0])

        # Initialize array sums to each corner
        array_sum = np.zeros((2, 2)).astype(float)
        for i in xrange(2):
            for j in xrange(2):
                removed_circle_array = remove_circle(image,
                                                     x[i+1], y[j+1],
                                                     radius, res=1)
                array_sum[j, i] = sum_array(removed_circle_array)

        # Find the index of the corner with minimum array_sum
        min_index = np.unravel_index(np.argmin(array_sum), (2, 2)) # Tuple

        while abs(x[3] - x[0]) > center_tol and abs(y[3] - y[0]) > center_tol:
            # Move the other corners to smaller search area
            if min_index[0] == 0:
                y[3] = y[2]
                y[2] = y[1]
                y[1] = y[0] + (1 - self._phi) * (y[3] - y[0])
            else:
                y[0] = y[1]
                y[1] = y[2]
                y[2] = y[0] + self._phi * (y[3] - y[0])
            if min_index[1] == 0:
                x[3] = x[2]
                x[2] = x[1]
                x[1] = x[0] + (1 - self._phi) * (x[3] - x[0])
            else:
                x[0] = x[1]
                x[1] = x[2]
                x[2] = x[0] + self._phi * (x[3] - x[0])

            # Replace the opposite corner array sum (so it doesn't need to be recalculated)
            array_sum[1 - min_index[0], 1 - min_index[1]] = array_sum[min_index]
            min_index = (1 - min_index[0], 1 - min_index[1])

            # Recalculate new sums for all four corners
            for i in xrange(2):
                for j in xrange(2):
                    if i != min_index[1] or j != min_index[0]:
                        temp_res = 1
                        if abs(x[3] - x[0]) < 10*center_tol and abs(y[3] - y[0]) < 10*center_tol:
                            temp_res = res
                        print x[i+1], y[j+1]
                        removed_circle_array = remove_circle(image,
                                                             x[i+1], y[j+1],
                                                             radius, temp_res)
                        array_sum[j, i] = sum_array(removed_circle_array)

            min_index = np.unravel_index(np.argmin(array_sum), (2, 2))

        self._center.circle.x = x[min_index[1]+1]
        self._center.circle.y = y[min_index[0]+1]
        self._diameter.circle = radius * 2.0
        self._array_sum.circle = np.amin(array_sum)

    def set_fiber_center_edge_method(self):
        """TAverages the fiber edges to set the fiber center

        Sets
        ----
        self._center.edge.y : float
        self._center.edge.x : float
        """
        self.set_fiber_edges()

        self._center.edge.y = (self._edges.top + self._edges.bottom) / 2.0
        self._center.edge.x = (self._edges.left + self._edges.right) / 2.0

    def set_fiber_edges(self):
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
        self._diameter.edge : float
        """
        filtered_image = self.get_filtered_image()

        left = -1
        right = -1
        for index in xrange(self.width):
            if left < 0:
                if filtered_image[:, index].max() > self.threshold:
                    left = index
            else:
                if filtered_image[:, index].max() > self.threshold:
                    right = index

        top = -1
        bottom = -1
        for index in xrange(self.height):
            if top < 0:
                if filtered_image[index, :].max() > self.threshold:
                    top = index
            else:
                if filtered_image[index, :].max() > self.threshold:
                    bottom = index

        self._edges.left = left
        self._edges.right = right
        self._edges.top = top
        self._edges.bottom = bottom
        self._diameter.edge = ((right - left) + (bottom - top)) / 2.0

    #=========================================================================#
    #==== Useful Methods =====================================================#
    #=========================================================================#

    def convert_fnum_to_radius(self, fnum, units):
        """Returns the value in the proper units"""
        return convert_fnum_to_radius(fnum,
                                      self.get_pixel_size(),
                                      self.get_magnification(),
                                      units)

    def plot_cross_sections(self, image_array=None, row=None, column=None):
        """Plots cross sections across the center of the fiber"""
        if image_array is None:
            image_array = self.get_image()
        if row is None:
            row = self.get_fiber_center()[0]
        if column is None:
            column = self.get_fiber_center()[1]
        plot_cross_sections(image_array, row, column)

#=============================================================================#
#===== Useful Functions ======================================================#
#=============================================================================#

def get_image_data(image_obj, **kwargs):
    """Returns relevant information from a FiberImage object

    Args
    ----
    image_obj : FiberImage
        Image object to be analyzed

    Returns
    -------
    image_array : ndarray
        the 2D image
    y0 : float
        the fiber center y
    x0 : float
        the fiber center x
    radius : float
        the fiber radius
    """
    y0, x0, diameter = image_obj.get_fiber_data(**kwargs)
    radius = diameter / 2.0
    image_array = image_obj.get_image()
    return image_array, y0, x0, radius

def convert_fnum_to_radius(fnum, pixel_size, magnification, units='pixels'):
    """Converts a focal ratio to an image radius in given units."""
    fcs_focal_length = 4.0 # inches
    diameter = fcs_focal_length / fnum # inches
    radius = 25400 * diameter / 2.0 # microns
    return convert_microns_to_units(radius, pixel_size, magnification, units)

if __name__ == "__main__":
    from .input_output import load_image_object

    folder = 'C:/Libraries/Box Sync/ExoLab/Fiber_Characterization/Image Analysis/data/scrambling/2016-08-05 Prototype Core Extension 1/'

    images = [folder + 'Shift_00/in_' + str(i).zfill(3) + '.fit' for i in xrange(10)]
    dark = [folder + 'Dark/in_' + str(i).zfill(3) + '.fit' for i in xrange(10)]
    ambient = [folder + 'Ambient/in_' + str(i).zfill(3) + '.fit' for i in xrange(10)]

    im_obj = FiberImage(images, dark=dark, ambient=ambient)

    tol = 0.25
    test_range = 5
    factor = 1.0

    im_obj.show_image_array()
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

    new_im_obj = load_image_object(im_obj.object_file, im_obj.image_file)

    new_im_obj.show_image_array()
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
