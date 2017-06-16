"""fiber_image.py was written by Ryan Petersburg for use with fiber
characterization on the EXtreme PREcision Spectrograph
"""
import numpy as np
from .numpy_array_handler import (sum_array, crop_image, remove_circle,
                                  isolate_circle, circle_array, polynomial_fit,
                                  gaussian_fit, rectangle_array,
                                  mesh_grid_from_array, intensity_array,
                                  isolate_rectangle)
from .plotting import (plot_cross_sections, plot_overlaid_cross_sections,
                       plot_dot, show_plots, plot_image)
from .containers import (FiberInfo, Edges, FRDInfo, ModalNoiseInfo,
                         convert_microns_to_units, Pixel)
from .calibrated_image import CalibratedImage
from .modal_noise import modal_noise

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
    def __init__(self, image_input, threshold=1000,
                 input_fnum=2.4, output_fnum=2.4, **kwargs):
        # Private attribute initialization
        self.threshold = threshold

        self._frd_info = FRDInfo()
        self._frd_info.input_fnum = input_fnum
        self._frd_info.output_fnum = output_fnum

        self._modal_noise_info = ModalNoiseInfo()

        self._phi = (5 ** 0.5 - 1) / 2

        self._gaussian_amp = 0.0
        self._gaussian_offset = 0.0

        super(FiberImage, self).__init__(image_input, **kwargs)

        self._edges = Edges(pixel_size=self.pixel_size,
                            magnification=self.magnification)
        self._center = FiberInfo('pixel', pixel_size=self.pixel_size,
                                 magnification=self.magnification)
        self._centroid = FiberInfo('pixel', pixel_size=self.pixel_size,
                                    magnification=self.magnification)
        self._diameter = FiberInfo('value')
        self._array_sum = FiberInfo('value')

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
        _center.method : fPixel
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

    def get_gaussian_fit(self, full_output=False, radius_factor=1.0):
        """Return the best gaussian fit for the image

        Returns
        -------
        _fit.gaussian : 2D numpy.ndarray
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

        if full_output:
            return gauss_fit, coeffs
        return gauss_fit

    def get_polynomial_fit(self, deg=6, radius_factor=0.95, fiber_method=None, **kwargs):
        """Return the best polynomial fit for the image

        Args
        ----
        deg : int (default=6)
            The degree of polynomial to fit

        Returns
        -------
        poly_fit : 2D numpy.ndarray
        """
        image = self.get_image()
        center = self.get_fiber_center(method=fiber_method, **kwargs)
        radius = self.get_fiber_radius() * radius_factor
        poly_fit = polynomial_fit(image, deg, center, radius)
        return poly_fit

    def get_tophat_fit(self):
        """Return the circle array that best covers the fiber face

        Returns
        -------
        circle_array : numpy.ndarray (2D)
            Circle array centered at best calculated center and with best
            calculated diameter
        """
        image = self.get_image()
        center = self.get_fiber_center()
        radius = self.get_fiber_radius()
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

    def set_fiber_centroid(self, method='full', radius_factor=1.0,
                           show_image=False, fiber_shape='circle', **kwargs):
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
        _centroid.method : Pixel
            The centroid of the image in the context of the given method
        """
        image = self.get_image()
        if method == 'full':
            image_iso = image * (self.get_filtered_image() > self.threshold).astype('float64')
        else:
            center = self.get_fiber_center(method=method, **kwargs)
            radius = self.get_fiber_radius(method=method, **kwargs)
            if 'rect' in fiber_shape:
                image_iso = isolate_rectangle(image, corners=self._edges)
            else:
                image_iso = isolate_circle(image, center,
                                           radius*radius_factor, res=1)

        x_array, y_array = self.get_mesh_grid()
        getattr(self._centroid, method).x = ((image_iso * x_array).sum()
                                             / image_iso.sum())
        getattr(self._centroid, method).y = ((image_iso * y_array).sum()
                                             / image_iso.sum())

        if show_image:
            plot_dot(image_iso, getattr(self._centroid, method))
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
        else:
            raise RuntimeError('Incorrect string for fiber centering method')

        if show_image:
            center = getattr(self._center, method)
            r = getattr(self._diameter, method) / 2.0
            image = self.get_filtered_image()

            if method == 'gaussian':
                plot_overlaid_cross_sections(image, self.get_gaussian_fit(),
                                             center)
                plot_dot(image, center)
                show_plots()
            else:
                plot_image(remove_circle(image, center, r, res=1))
                plot_overlaid_cross_sections(image, image.max() / 2.0
                                             * circle_array(self.get_mesh_grid(),
                                                            center.x, center.y, r, res=1),
                                             center)
                if method == 'edge':
                    for corner in self._edges:
                        plot_dot(image, corner)
                show_plots()

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
        _, coeffs = self.get_gaussian_fit(full_output=True)

        self._center.gaussian.x = coeffs[0]
        self._center.gaussian.y = coeffs[1]
        self._diameter.gaussian = abs(coeffs[2]) * 2.0
        self._gaussian_amp = coeffs[3]
        self._gaussian_offset = coeffs[4]

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
            approx_radius = self.get_fiber_radius(method='edge')
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
                                       center_range=None, image=None, **kwargs):
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
            radius = self.get_fiber_radius(method='edge')

        # Create four "corners" to test center of the removed circle
        x = np.zeros(4).astype(float)
        y = np.zeros(4).astype(float)

        if center_range is not None:
            approx_center = self.get_fiber_center(method='edge')
            center_range = center_range / 2.0

            x[0] = approx_center.x - center_range
            if x[0] < radius:
                x[0] = radius
            x[3] = approx_center.x + center_range
            if x[3] > self.width - radius:
                x[3] = self.width - radius

            y[0] = approx_center.y - center_range
            if y[0] < radius:
                y[0] = radius
            y[3] = approx_center.y + center_range
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
                                                     Pixel(x[i+1], y[j+1]),
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
                        removed_circle_array = remove_circle(image,
                                                             Pixel(x[i+1], y[j+1]),
                                                             radius, temp_res)
                        array_sum[j, i] = sum_array(removed_circle_array)

            min_index = np.unravel_index(np.argmin(array_sum), (2, 2))

        self._center.circle.x = x[min_index[1]+1]
        self._center.circle.y = y[min_index[0]+1]
        self._diameter.circle = radius * 2.0
        self._array_sum.circle = np.amin(array_sum)

    def set_fiber_center_edge_method(self, **kwargs):
        """TAverages the fiber edges to set the fiber center

        Sets
        ----
        self._center.edge.y : float
        self._center.edge.x : float
        """
        self.set_fiber_edges(**kwargs)

        self._center.edge.y = (self._edges.top.y + self._edges.bottom.y) / 2.0
        self._center.edge.x = (self._edges.left.x + self._edges.right.x) / 2.0

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
        self._diameter.edge : float
        """
        image = self.get_filtered_image() # To prvent hot pixels

        left = -1
        right = -1
        for index in xrange(self.width):
            if left < 0:
                if image[:, index].max() > self.threshold:
                    left = index
            else:
                if image[:, index].max() > self.threshold:
                    right = index

        top = -1
        bottom = -1
        for index in xrange(self.height):
            if top < 0:
                if image[index, :].max() > self.threshold:
                    top = index
            else:
                if image[index, :].max() > self.threshold:
                    bottom = index

        left = np.array([left, image[:, left].argmax()])
        right = np.array([right, image[:, right].argmax()])
        top = np.array([image[top, :].argmax(), top])
        bottom = np.array([image[bottom, :].argmax(), bottom])
        diameter = (np.sqrt(((right - left)**2).sum())
                  + np.sqrt(((bottom - top)**2).sum())) / 2.0

        self._edges.left.set_pixel(*left)
        self._edges.right.set_pixel(*right)
        self._edges.top.set_pixel(*top)
        self._edges.bottom.set_pixel(*bottom)
        self._diameter.edge = diameter

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
