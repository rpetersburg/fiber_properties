"""scrambling_gain.py was written by Ryan Petersburg for use with fiber
characterization for the EXtreme PRecision Spectrograph

This module contains functions that calculate the scrambling gain for
multiple FCS images contained in FiberImage objects
"""
from collections import Iterable
import numpy as np
from .fiber_image import FiberImage
from .containers import ScramblingInfo

def scrambling_gain(in_objs, out_objs, input_method=None, output_method=None, **kwargs):
    """Calculates the scrambling gain for fiber input and output images

    Args
    ----
    in_objs : list(FiberImage) or list(str)
        list of the FiberImage input objects or saved object file names
    out_objs : list(FiberImage) or list(str)
        list of the FiberImage output objects or saved object file names
    input_method : str {'edge', 'radius'}, optional
        method used to find the diameter of the input fiber face
    output_method : str {'edge','radius','gaussian'}, optional
        method used to find the diameter of the output fiber image

    Returns
    -------
    info : ScramblingInfo()
        Object containing scrambling information. See containers.py for
        specifics
    """
    if not isinstance(in_objs, Iterable):
        in_objs = [in_objs]
    if not isinstance(out_objs, Iterable):
        out_objs = [out_objs]

    if len(in_objs) != len(out_objs):
        raise RuntimeError('Lists of input and output objects not the same length')

    info = ScramblingInfo()

    for in_obj in in_objs:
        if isinstance(in_obj, basestring):
            in_obj = FiberImage(in_obj)
        in_centroid = in_obj.get_fiber_centroid(method=input_method,
                                                units='microns',
                                                **kwargs)
        in_center = in_obj.get_fiber_center(method=input_method,
                                            units='microns',
                                            **kwargs)
        in_diameter = in_obj.get_fiber_diameter(method=input_method,
                                                units='microns',
                                                **kwargs)
        in_obj.save()
        info.in_x.append((in_centroid.x - in_center.x) / in_diameter)
        info.in_y.append((in_centroid.y - in_center.y) / in_diameter)

    for out_obj in out_objs:
        if isinstance(out_obj, basestring):
            out_obj = FiberImage(out_obj)
        out_centroid = out_obj.get_fiber_centroid(method=output_method,
                                                  units='microns',
                                                  **kwargs)
        out_center = out_obj.get_fiber_center(method=output_method,
                                              units='microns',
                                              **kwargs)
        out_diameter = out_obj.get_fiber_diameter(method=output_method,
                                                  units='microns',
                                                  **kwargs)
        out_obj.save()
        info.out_x.append((out_centroid.x - out_center.x) / out_diameter)
        info.out_y.append((out_centroid.y - out_center.y) / out_diameter)

    list_len = len(info.in_x)
    for i in xrange(list_len):
        for j in xrange(i+1, list_len):
            info.in_d.append(np.sqrt((info.in_x[i] - info.in_x[j])**2
                                     + (info.in_y[i] - info.in_y[j])**2))
            info.out_d.append(np.sqrt((info.out_x[i] - info.out_x[j])**2
                                      + (info.out_y[i] - info.out_y[j])**2))
    info.scrambling_gain = np.array(info.in_d) / np.array(info.out_d)
    # for i in xrange(list_len):
    #     d_in = np.sqrt((info.in_x[i] - info.in_x[0])**2 + (info.in_y[i] - info.in_y[0])**2)
    #     d_out = np.sqrt((info.out_x[i] - info.out_x[0])**2 + (info.in_y[i] - info.in_y[0])**2)
    #     scrambling_gain.append(d_in / d_out)

    return info
