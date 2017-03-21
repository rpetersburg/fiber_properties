"""scrambling_gain.py was written by Ryan Petersburg for use with fiber
characterization for the EXtreme PRecision Spectrograph

This module contains functions that calculate the scrambling gain for
multiple FCS images contained in ImageAnalysis objects
"""
from collections import Iterable
import numpy as np
from fiber_properties.input_output import load_image_object
from fiber_properties.containers import ScramblingInfo

def scrambling_gain(in_objs, out_objs, input_method=None, output_method=None, **kwargs):
    """Calculates the scrambling gain for fiber input and output images

    Args
    ----
    in_objs : list(ImageAnalysis) or list(string)
        list of the ImageAnalysis input objects or object file names
    out_objs : list(ImageAnalysis) or list(string)
        list of the ImageAnalysis output objects or object file names
    input_method : str {'edge', 'radius'}, optional
        method used to find the diameter of the input fiber face
    output_method : str {'edge','radius','gaussian'}, optional
        method used to find the diameter of the output fiber image

    Returns
    -------
    info : ScramblingInfo()
        Object containing scrambling information. See
        fiber_properties.containers for specifics
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
            in_obj = load_image_object(in_obj)
        in_centroid = in_obj.get_fiber_centroid(radius_factor=1.05,
                                                method='gaussian',
                                                units='microns',
                                                **kwargs)
        in_center = in_obj.get_fiber_center(method=input_method,
                                            units='microns',
                                            **kwargs)
        in_diameter = in_obj.get_fiber_diameter(method=input_method,
                                                units='microns',
                                                **kwargs)
        in_obj.save()
        info.in_x.append((in_centroid[1] - in_center[1]) / in_diameter)
        info.in_y.append((in_centroid[0] - in_center[0]) / in_diameter)

    for out_obj in out_objs:
        if isinstance(out_obj, basestring):
            out_obj = load_image_object(out_obj)
        out_centroid = out_obj.get_fiber_centroid(radius_factor=1.0,
                                                  method=output_method,
                                                  units='microns',
                                                  **kwargs)
        out_center = out_obj.get_fiber_center(method=output_method,
                                              units='microns',
                                              **kwargs)
        out_diameter = out_obj.get_fiber_diameter(method=output_method,
                                                  units='microns',
                                                  **kwargs)
        out_obj.save()
        info.out_x.append((out_centroid[1] - out_center[1]) / out_diameter)
        info.out_y.append((out_centroid[0] - out_center[0]) / out_diameter)

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
