"""scrambling_gain.py was written by Ryan Petersburg for use with fiber
characterization for the EXtreme PRecision Spectrograph

This module contains functions that calculate the scrambling gain for
multiple FCS images contained in ImageAnalysis objects
"""
import numpy as np
from collections import Iterable
from input_output import load_image_object

def scrambling_gain(in_objs, out_objs, input_method=None, output_method=None):
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
    input_x : list
    input_y : list
    output_x : list
    output_y : list        
    scrambling_gain : list
    input_dist : list
    output_dist : list
    """
    if not isinstance(in_objs, Iterable):
        in_objs = [in_objs]
    if not isinstance(out_objs, Iterable):
        out_objs = [out_objs]

    if len(in_objs) != len(out_objs):
        raise RuntimeError('Lists of input and output objects not the same length')

    input_x = []
    input_y = []
    for in_obj in in_objs:
        if isinstance(in_obj, basestring):
            in_obj = load_image_object(in_obj)
        in_centroid = in_obj.getFiberCentroid(radius_factor=1.05, method='gaussian', units='microns')
        in_center = in_obj.getFiberCenter(method=input_method, units='microns')
        in_diameter = in_obj.getFiberDiameter(method=input_method, units='microns')
        in_obj.save()
        input_x.append((in_centroid[1] - in_center[1]) / in_diameter)
        input_y.append((in_centroid[0] - in_center[0]) / in_diameter)

    output_x = []
    output_y = []
    for out_obj in out_objs:
        if isinstance(out_obj, basestring):
            out_obj = load_image_object(out_obj)
        out_centroid = out_obj.getFiberCentroid(radius_factor=1.0, method=output_method, units='microns')
        out_center = out_obj.getFiberCenter(method=output_method, units='microns')
        out_diameter = out_obj.getFiberDiameter(method=output_method, units='microns')
        out_obj.save()
        output_x.append((out_centroid[1] - out_center[1]) / out_diameter)
        output_y.append((out_centroid[0] - out_center[0]) / out_diameter)

    scrambling_gain = []
    input_dist = []
    output_dist = []
    list_len = len(input_x)
    for i in xrange(list_len):
        for j in xrange(i+1, list_len):
            d_in.append(np.sqrt((input_x[i] - input_x[j])**2 + (input_y[i] - input_y[j])**2))
            d_out.append(np.sqrt((output_x[i] - output_x[j])**2 + (output_y[i] - output_y[j])**2))
            scrambling_gain.append(d_in[i+j] / d_out[i+j])
    # for i in xrange(list_len):
    #     d_in = np.sqrt((input_x[i] - input_x[0])**2 + (input_y[i] - input_y[0])**2)
    #     d_out = np.sqrt((output_x[i] - output_x[0])**2 + (input_y[i] - input_y[0])**2)
    #     scrambling_gain.append(d_in / d_out)

    return input_x, input_y, output_x, output_y, scrambling_gain, d_in, d_out