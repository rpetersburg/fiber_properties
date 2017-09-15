"""fRatioDegradation.py was written by Ryan Petersburg for use with fiber
characterization for the EXtreme PRecision Spectrograph

This module contains functions that calculate the focal ratio degradation
for images taken with the FCS
"""
import numpy as np
import scipy.stats as stats
from .containers import FRDInfo
from .input_output import load_image_object
from .fiber_image import FiberImage

def frd(in_objs, out_objs, cal_method='full', save_objs=True, **kwargs):
    """Collects all relevant FRD info from the frd_input

    Args
    ----
    in_objs : list(FiberImage) or list(str)
        List of objects (or file names of saved objects) containing far field
        images taken at different input focal ratios
    out_objs : list(FiberImage) or list(str)
        List of objects (or file names of saved objects) containing far field
        images taken at different output focal ratios
    cal_method : str, optional
        Method used to calculate the diameter of the output images
    save_objs : bool, optional
        If true, the FiberImage objects will be saved after calculations
        are made
    **kwargs : **dict
        Keyword arguments that are passed to FiberImage.get_frd_info

    Returns
    -------
    output : FRDInfo
        object containing frd information. See fiber_properties.containers
        for specifics
    magnification : float
        the averaged magnification for the far field camera
    magn_list : list(float)
        a list of each measured magnification for the far field camera
    magn_error : float
        the relative error in the focal ratios as calculated from the
        magnification values
    """
    output = FRDInfo()

    magn_list = []
    for out_obj in out_objs:
        if isinstance(out_obj, basestring):
            out_obj = FiberImage(out_obj)
        diameter = out_obj.get_fiber_diameter(method=cal_method,
                                              units='microns')
        if save_objs:
            out_obj.save_object()
        magn_list.append(diameter / ((4.0 / out_obj.get_output_fnum()) * 25400))

    magnification = np.mean(magn_list)
    magn_error = 0.0
    if len(magn_list) > 1:
        magn_error = stats.sem(magn_list)

    for in_obj in in_objs:
        if isinstance(in_obj, basestring):
            in_obj = FiberImage(in_obj)
        in_obj.set_magnification(magnification)

        temp_output = in_obj.get_frd_info(**kwargs)

        if save_objs:
            in_obj.save_object()

        for attr in vars(temp_output):
            getattr(output, attr).append(getattr(temp_output, attr))

    return output, magnification, magn_list, magn_error
