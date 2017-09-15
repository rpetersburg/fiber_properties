"""fiber_centroid.py was written by Ryan Petersburg for use with fiber
characterization on the EXtreme PREcision Spectrograph
The functions in this module calculate the centroid of a FiberImage object
"""
from .plotting import plot_dot, show_plots
from .containers import Pixel
from .numpy_array_handler import isolate_rectangle, isolate_circle, mesh_grid_from_array

def fiber_centroid(im_obj, method='full', radius_factor=1.0,
                   show_image=False, fiber_shape='circle', **kwargs):
    """Find the centroid of the fiber face image
    Args
    ----
    method : {'full', 'edge', 'radius', 'gaussian', 'circle'}, optional
        If 'full', takes the centroid of the entire image. Otherwise, uses
        the specified method to isolate only the fiber face in the image
    radius_factor : number
        The factor by which the radius is multiplied when isolating the
        fiber face in the image
    show_image : bool, optional
        Shows centroid dot on the fiber image
    fiber_shape : {'circle', 'rectangle'}, optional
        The shape of the fiber core cross-section. Used to decide which
        points to use when calculating the centroid.
    **kwargs :
        The keyworded arguments to pass to the centering method
    Returns
    -------
    centroid : Pixel
        Centroid of the fiber face image
    """
    image = im_obj.get_image()
    if method == 'full':
        image_iso = image
    else:
        center = im_obj.get_fiber_center(method=method, **kwargs)
        radius = im_obj.get_fiber_radius(method=method, **kwargs)
        if 'rect' in fiber_shape:
            image_iso = isolate_rectangle(image, corners=self._edges)
        else:
            image_iso = isolate_circle(image, center,
                                       radius*radius_factor, res=1)
    image_iso *= (im_obj.get_filtered_image()
                  > im_obj.threshold).astype('uint8')

    centroid = calc_centroid(image_iso)

    if show_image:
        plot_dot(image, centroid.as_tuple()[::-1])
        show_plots()

    return centroid

def calc_centroid(image):
    """Calculates the 2D centroid of a full image."""
    x_array, y_array = mesh_grid_from_array(image)
    centroid = Pixel(units='pixels')
    centroid.x = (image * x_array).sum() / image.sum()
    centroid.y = (image * y_array).sum() / image.sum()
    return centroid