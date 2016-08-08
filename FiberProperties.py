import numpy as np
import matplotlib.pyplot as plt

#=============================================================================#
#==== Scrambling Gain Functions ==============================================#
#=============================================================================#

def scramblingGain(in_obj_1, out_obj_1, in_obj_2=None, out_obj_2=None, method='edge'):
    """Calculates the scrambling gain for fiber input and output images

    Args:
        in_obj_1: The first input object


    """
    in_centroid_1 = in_obj_1.getFiberCentroid()
    in_diameter_1 = in_obj_1.getFiberDiameter(method=method)

    out_centroid_1 = out_obj_1.getFiberCentroid()
    out_diameter_1 = out_obj_1.getFiberDiameter(method=method)

    if in_obj_2 is None or out_obj_2 is None:
        in_centroid_2 = in_obj_1.getFiberCenter(method=method)
        out_centroid_2 = out_obj_1.getFiberCenter(method=method)
        in_diameter = in_diameter_1
        out_diameter = out_diameter_1

    else:
        in_centroid_2 = in_obj_2.getFiberCentroid()
        out_centroid_2 = out_obj_2.getFiberCentroid()
        in_diameter = (in_diameter_1 + in_obj_2.getFiberDiameter(method=method)) / 2.0
        out_diameter = (out_diameter_1 + in_obj_2.getFiberDiameter(method=method)) / 2.0


    delta_D_in = np.sqrt((in_centroid_1[1] - in_centroid_2[1])**2 + (in_centroid_1[0] - in_centroid_2[0])**2)
    delta_D_out = np.sqrt((out_centroid_1[1] - out_centroid_2[1])**2 + (out_centroid_1[0] - out_centroid_2[0])**2)

    scramblingGain = (delta_D_in / in_diameter) / (delta_D_out / out_diameter)

    return scramblingGain

#=============================================================================#
#==== Focal Ratio Degradation Functions ======================================#
#=============================================================================#

def FRD(ff_obj, input_f_number=-1.0, f_lim=(2.4, 10.0), res=0.1):
    center_y, center_x = ff_obj.getFiberCentroid()

    f_numbers = np.arange(f_lim[0], f_lim[1], res)
    encircled_energy = []
    for f in f_numbers:
        radius = fnumber_to_radius(f, ff_obj.magnification)
        isolated_circle = ff_obj.isolateCircle(ff_obj.getImage(),
                                               center_x,
                                               center_y,
                                               radius)
        iso_circ_sum = isolated_circle.sum()
        encircled_energy.append(iso_circ_sum)
        if abs(f - input_f_number) < res / 2.0:
            energy_loss = 100 * (1 - iso_circ_sum / encircled_energy[0])
        if iso_circ_sum / encircled_energy[0] >= 0.95:
            output_f_number = f

    encircled_energy = list(np.array(encircled_energy) / encircled_energy[0])

    return f_numbers, encircled_energy, energy_loss, output_f_number


def fnumber_to_radius(fnumber, magnification):
    return magnification * 4.0 * 25400 / fnumber / 2.0 / 3.45
