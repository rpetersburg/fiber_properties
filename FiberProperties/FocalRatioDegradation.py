"""FocalRatioDegradation.py was written by Ryan Petersburg for use with fiber
characterization for the EXtreme PRecision Spectrograph

This module contains functions that calculate the focal ratio degradation
for images taken with the FCS contained in ImageAnalysis objects
"""

from NumpyArrayHandler import isolateCircle

def FRD(ff_obj, input_focal_ratio=-1.0, focal_lim=(2.4, 10.0), res=0.1):
    """Calculates the encircled energy for various f ratios

    Args:
        ff_obj [ImageAnalysis]: the image object on which the FRD is calculated
        input_focal_ratio [float]: the fiber input f ratio
        focal_lim [(float, float)]: the limits of the f ratio calculations
        res [float]: the spacing between each f ratio when calculating
            encircld energy

    Returns:
        focal_ratios [list]: list of the f ratios used to calculate encircled
            energy
        encircled_energy [list]: list of the encircled energy at each given
            f ratio
        energy_loss [float]: the loss of energy when the output f ratio equals
            the input f ratio given as a percent
        output_focal_ratio [float]: the approximate f ratio inside which 95% of 
            the total encircled energy is included
    """
    center_y, center_x = ff_obj.getFiberCentroid()

    focal_ratios = list(np.arange(focal_lim[0], focal_lim[1] + res, res))
    energy_loss = None
    output_focal_ratio = None
    encircled_energy = []
    for f in focal_ratios:
        radius = _focal_ratio_to_radius(f, ff_obj)
        isolated_circle = isolateCircle(ff_obj.getImage(),
                                        center_x,
                                        center_y,
                                        radius)
        iso_circ_sum = isolated_circle.sum()
        encircled_energy.append(iso_circ_sum)
        if abs(f - input_focal_ratio) < res / 2.0:
            energy_loss = 100 * (1 - iso_circ_sum / encircled_energy[0])
        if iso_circ_sum / encircled_energy[0] >= 0.95:
            output_focal_ratio = f

    encircled_energy = list(np.array(encircled_energy) / encircled_energy[0])

    return focal_ratios, encircled_energy, energy_loss, output_focal_ratio

def _focal_ratio_to_radius(focal_ratio, im_obj):
    """Converts an f ratio to an image radius in units of pixels
    
    Args:
        focal_ratio [float]: the f ratio to be converted
        magnification [float]: the magnification of the camera

    Returns:
        radius [float]: in units of pixels
    """
    return 25400 * (4.0 / focal_ratio) * (im_obj.magnification / im_obj.pixel_size) / 2.0