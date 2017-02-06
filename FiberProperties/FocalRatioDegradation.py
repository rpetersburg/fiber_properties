"""FocalRatioDegradation.py was written by Ryan Petersburg for use with fiber
characterization for the EXtreme PRecision Spectrograph

This module contains functions that calculate the focal ratio degradation
for images taken with the FCS contained in ImageAnalysis objects
"""
import cPickle as pickle
import numpy as np
from scipy import stats
from NumpyArrayHandler import isolateCircle
from Calibration import Calibration
from ImageAnalysis import ImageAnalysis

FRD_CALIBRATION_THRESHOLD = 1500
FOCAL_RATIO_DIAMETER = 0.95

class FRD_Input(object):
    """Container for FRD test input information

    Attributes
    ----------
    name : string
        Name to be used for saving plots
    folder : string
        Top level folder where images are contained
    input_focal_ratios : list(float), optional
        Input focal ratios which have associated images
    cal_focal_ratios : list(float), optional
        Output focal ratios which were used as calibration images
    """
    def __init__(self, name, folder,
                 input_focal_ratios=[2.5, 3.0, 3.5, 4.0, 4.5, 5.0],
                 cal_focal_ratios=[3.0, 4.0, 5.0]):
        self.name = name
        self.folder = folder
        self.input_focal_ratios = input_focal_ratios
        self.cal_focal_ratios = cal_focal_ratios

class FRD_Output(object):
    """Container for FRD test output information

    Attributes
    ----------
    input_focal_ratio : list(float)
        list of the given input focal ratios
    encircled_energy : list(list(float))
        list of the lists of encircled energies for each input focal
        ratio
    encircled_energy_focal_ratios : list(list(float))
        independent variable (output f/#) corresponding to each
        encircled energy
    energy_loss : list(float)
        list of energy losses for each input focal ratio
    output_focal_ratios : list(float)
        list of calculated output focal ratio for each input focal
        ratio
    magnification : float
                the averaged magnification for the far field camera
    magn_list : list(float)
        a list of each measured magnification for the far field camera
    magn_error : float
        the relative error in the focal ratios as calculated from the
        magnification values
    """
    def __init__(self):
        self.encircled_energy_focal_ratios = []
        self.encircled_energy = []
        self.output_focal_ratios = []
        self.energy_loss = []
        self.magnification = None
        self.magn_list = []
        self.magn_error = 0.0

def FRD(frd_input, save_images=True):
    """Collects all relevant FRD info from the frd_input

    Args
    ----
    frd_input : FRD_Input
        object containing test information

    Returns
    -------
    frd_output : FRD_Output
        object containing frd information
    """
    print 'Starting ' + frd_input.name
    folder = frd_input.folder
    input_focal_ratios = frd_input.input_focal_ratios
    cal_focal_ratios = frd_input.cal_focal_ratios
    frd_output = FRD_Output()

    calibration = Calibration(dark=[folder + 'Dark/im_' + str(i).zfill(3) + '.fit' for i in xrange(10)],
                              ambient=[folder + 'Ambient/im_' + str(i).zfill(3) + '.fit' for i in xrange(10)])

    for f in cal_focal_ratios:
        images = [folder + 'Output ' + str(f) + '/im_' + str(i).zfill(3) + '.fit' for i in xrange(10)]
        im_obj = ImageAnalysis(images, calibration, magnification=1,
                               threshold=FRD_CALIBRATION_THRESHOLD)
        diameter = im_obj.getFiberDiameter(method='edge',
                                           units='microns',
                                           show_image=False)
        frd_output.magn_list.append(diameter / ((4.0 / f) * 25400))

        if save_images:
            im_obj.saveImage(folder + 'Output ' + str(f) + ' Image.fit')

    frd_output.magnification = np.array(frd_output.magn_list).mean()
    if len(frd_output.magn_list) > 1:
        frd_output.magn_error = stats.sem(frd_output.magn_list)

    print frd_input.name + ' calibration complete'

    for f in input_focal_ratios:
        images = [folder + 'Input ' + str(f) + '/im_' + str(i).zfill(3) + '.fit' for i in xrange(10)]
        im_obj = ImageAnalysis(images, calibration,
                               magnification=frd_output.magnification)
        temp_output = singleImageFRD(im_obj, input_focal_ratio=f,
                                     focal_lim=(2.3, 6.0), res=0.1)

        frd_output.encircled_energy_focal_ratios.append(temp_output[0])
        frd_output.encircled_energy.append(temp_output[1])
        frd_output.energy_loss.append(temp_output[2])
        frd_output.output_focal_ratios.append(temp_output[3])

        if save_images:
            im_obj.saveImage(folder + 'Input ' + str(f) + ' Image.fit')

    with open(folder + 'FRD_Output.pkl', 'wb') as output_file:
        pickle.dump(frd_output, output_file)
    with open(folder + 'FRD_Output.txt', 'w') as output_file:
        output_file.write(str(frd_output.__dict__))

    print frd_input.name + ' FRD calculations complete'

    return frd_output

def singleImageFRD(im_obj, input_focal_ratio=-1.0,
                   focal_lim=(2.3, 6.0), res=0.1):
    """Calculate the encircled energy for various f ratios

    Args
    ----
    im_obj : ImageAnalysis
        the image object on which the FRD is calculated
    input_focal_ratio :float
        the fiber input f ratio
    focal_lim : (float, float)
        the limits of the f ratio calculations
    res : float
        the spacing between each f ratio when calculating encircled energy

    Returns
    -------
    focal_ratios : list(float)
        list of the f ratios used to calculate encircled energy
    encircled_energy : list(float)
        list of the encircled energy at each given f ratio
    energy_loss : float
        the loss of energy when the output f ratio equals the input f ratio
        given as a percent
    output_focal_ratio : float
        the approximate f ratio inside which 95% of the total encircled energy
        is included
    """
    center_y, center_x = im_obj.getFiberCentroid()

    focal_ratios = list(np.arange(focal_lim[0], focal_lim[1] + res, res))
    energy_loss = None
    output_focal_ratio = None
    encircled_energy = []
    for f in focal_ratios:
        radius = _focal_ratio_to_radius(f, im_obj)
        isolated_circle = isolateCircle(im_obj.getImage(),
                                        center_x,
                                        center_y,
                                        radius)
        iso_circ_sum = isolated_circle.sum()
        encircled_energy.append(iso_circ_sum)
        if abs(f - input_focal_ratio) < res / 2.0:
            energy_loss = 100 * (1 - iso_circ_sum / encircled_energy[0])
        if iso_circ_sum / encircled_energy[0] >= FOCAL_RATIO_DIAMETER:
            output_focal_ratio = f

    encircled_energy = list(np.array(encircled_energy) / encircled_energy[0])

    return focal_ratios, encircled_energy, energy_loss, output_focal_ratio

def _focal_ratio_to_radius(focal_ratio, im_obj):
    """Converts an f ratio to an image radius in units of pixels

    Args
    ----
    focal_ratio : float
        the f ratio to be converted
    magnification : float
        the magnification of the camera

    Returns
    -------
    radius : float
        in units of pixels
    """
    return 25400 * (4.0 / focal_ratio) * (im_obj.getMagnification() / im_obj.getPixelSize()) / 2.0
