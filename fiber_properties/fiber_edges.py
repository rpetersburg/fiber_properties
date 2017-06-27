"""fiber_edges.py was written by Ryan Petersburg for use with fiber
characterization on the EXtreme PREcision Spectrograph

The functions in this module are used to calculate the location of the edges of
a fiber on a FiberImage object
"""

import numpy as np
from .containers import Edges

def fiber_edges(im_obj, **kwargs):
    """Return the fiber edge pixel values

    Finds the left, right, top, and bottom edges of the fiber by finding where
    the maxima of each row and column cross the given threshold. Also sets
    the width of the fiber by the maximum of the horizontal and vertical
    lengths

    Returns
    -------
    edges : .containers.Edges
        container with the edge Pixels
    """
    image = im_obj.get_filtered_image() # To prvent hot pixels

    left = -1
    right = -1
    for index in xrange(im_obj.width):
        if left < 0:
            if image[:, index].max() > im_obj.threshold:
                left = index
        else:
            if image[:, index].max() > im_obj.threshold:
                right = index

    top = -1
    bottom = -1
    for index in xrange(im_obj.height):
        if top < 0:
            if image[index, :].max() > im_obj.threshold:
                top = index
        else:
            if image[index, :].max() > im_obj.threshold:
                bottom = index

    left = np.array([left, image[:, left].argmax()])
    right = np.array([right, image[:, right].argmax()])
    top = np.array([image[top, :].argmax(), top])
    bottom = np.array([image[bottom, :].argmax(), bottom])

    edges = Edges()
    edges.left.set_pixel(*left)
    edges.right.set_pixel(*right)
    edges.top.set_pixel(*top)
    edges.bottom.set_pixel(*bottom)

    return edges
