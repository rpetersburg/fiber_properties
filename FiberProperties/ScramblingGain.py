"""ScramblingGain.py was written by Ryan Petersburg for use with fiber
characterization for the EXtreme PRecision Spectrograph

This module contains functions that calculate the scrambling gain for
multiple FCS images contain in ImageAnalysis objects
"""
import numpy as np

def scramblingGain(in_objs, out_objs, input_method=None, output_method=None):
    """Calculates the scrambling gain for fiber input and output images

    Args:
        in_objs : list(ImageAnalysis)
            list of the ImageAnalysis input objects
        out_objs : list(ImageAnalysis)
            list of the ImageAnalysis output objects
        input_method : str {'edge', 'radius'}, optional
            method used to find the diameter of the input fiber face
        output_method : str {'edge','radius','gaussian'}, optional
            method used to find the diameter of the output fiber image

    Returns:
        input_x [list]
        input_y [list]
        output_x [list]
        output_y [list]        
        scrambling_gain [float]: approximate scrambling gain
    """
    if not isinstance(in_objs, Iterable) or not isinstance(out_objs, Iterable):
        in_objs = [in_objs]
        out_objs = [out_objs]

    if len(in_objs) != len(out_objs):
        raise RuntimeError('Lists of input and output objects not the same length')

    input_x = []
    input_y = []
    for in_obj in in_objs:
        in_centroid = in_obj.getFiberCentroid(radius_factor=1.05, method='gaussian', units='microns')
        in_center = in_obj.getFiberCenter(method=input_method, units='microns')
        in_diameter = in_obj.getFiberDiameter(method=input_method, units='microns')
        input_x.append((in_centroid[1] - in_center[1]) / in_diameter)
        input_y.append((in_centroid[0] - in_center[0]) / in_diameter)

    output_x = []
    output_y = []
    for out_obj in out_objs:
        out_centroid = out_obj.getFiberCentroid(radius_factor=1.0, method=output_method, units='microns')
        out_center = out_obj.getFiberCenter(method=output_method, units='microns')
        out_diameter = out_obj.getFiberDiameter(method=output_method, units='microns')
        output_x.append((out_centroid[1] - out_center[1]) / out_diameter)
        output_y.append((out_centroid[0] - out_center[0]) / out_diameter)

    scrambling_gain = []
    list_len = len(input_x)
    for i in xrange(list_len):
        for j in xrange(i+1, list_len):
            d_in = np.sqrt((input_x[i] - input_x[j])**2 + (input_y[i] - input_y[j])**2)
            d_out = np.sqrt((output_x[i] - output_x[j])**2 + (output_y[i] - output_y[j])**2)
            scrambling_gain.append(d_in / d_out)

    return input_x, input_y, output_x, output_y, scrambling_gain