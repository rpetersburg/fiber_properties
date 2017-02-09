"""fRatioDegradation.py was written by Ryan Petersburg for use with fiber
characterization for the EXtreme PRecision Spectrograph

This module contains functions that calculate the focal ratio degradation
for images taken with the FCS
"""
from Containers import FRDInfo
from InputOutput import loadImageObject

def FRD(in_objs, out_objs, cal_method='edge', save_objs=True, **kwargs):
    """Collects all relevant FRD info from the frd_input

    Args
    ----
    frd_input : FRD_Input
        object containing test information

    Returns
    -------
    output : FRD_Output
        object containing frd information
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
            out_obj = loadImageObject(out_obj)
        diameter = out_obj.getFiberDiameter(method=cal_method,
                                            units='microns')
        if save_objs:
            out_obj.saveObject()
        magn_list.append(diameter / ((4.0 / out_obj.getOutputFnum()) * 25400))

    magnification = np.array(output.magn_list).mean()
    if len(magn_list) > 1:
        magn_error = stats.sem(output.magn_list)

    for in_obj in in_objs:
        if isinstance(in_obj, basestring):
            in_obj = loadImageObject(in_obj)
        in_obj.setMagnification(magnification)
        
        temp_output = in_obj.getFRDInfo(**kwargs)
        if save_objs:
            in_obj.saveObject()

        for attr in vars(temp_output):
            getattr(output, attr).append(getattr(temp_output, attr))

    return output, magnification, magn_list, magn_error
