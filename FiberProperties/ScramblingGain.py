"""ScramblingGain.py was written by Ryan Petersburg for use with fiber
characterization for the EXtreme PRecision Spectrograph

This module contains functions that calculate the scrambling gain for
multiple FCS images contain in ImageAnalysis objects
"""

def scramblingGain(in_objs, out_objs, in_method='gaussian', out_method='edge'):
    """Calculates the scrambling gain for fiber input and output images

    Args:
        in_objs [list(ImageAnalysis)]: list of the ImageAnalysis input objects
        out_objs [list(ImageAnalysis)]: list of the ImageAnalysis output objects

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
        in_centroid = in_obj.getFiberCentroid(radius_factor=1.05, method=in_method, units='microns')
        in_center = in_obj.getFiberCenter(method=out_method, units='microns')
        in_diameter = in_obj.getFiberDiameter(method=out_method, units='microns')
        input_x.append((in_centroid[1] - in_center[1]) / in_diameter)
        input_y.append((in_centroid[0] - in_center[0]) / in_diameter)

    output_x = []
    output_y = []
    output_diameter = out_objs[0].getFiberDiameter(method=out_method, units='microns')
    for out_obj in out_objs:
        out_centroid = out_obj.getFiberCentroid(radius_factor=1.0, method=out_method, units='microns')
        out_center = out_obj.getFiberCenter(method=out_method, units='microns')
        out_diameter = out_obj.getFiberDiameter(method=out_method, units='microns')
        output_x.append((out_centroid[1] - out_center[1]) / out_diameter)
        output_y.append((out_centroid[0] - out_center[0]) / out_diameter)

    scrambling_gain = None

    return input_x, input_y, output_x, output_y, scrambling_gain

def twoPointScramblingGain(in_centroid_1, in_centroid_2, out_centroid_1, out_centroid_2):
    return