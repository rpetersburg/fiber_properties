from FiberProperties import ImageAnalysis, Calibration, FRD, saveArray, saveCrossSections
import matplotlib.pyplot as plt
import numpy as np
from ast import literal_eval
from sys import platform

class Data(object):
    """Container for relevant test information

    Attributes
    ----------
    name : string
        Name to be used for saving plots
    folder : string
        Top level folder where images are contained
    new : boolean
        If True, obtain new data from raw images. If False, open the saved
        data file in the top level directory
    input_focal_ratios : list(float), optional
        Input focal ratios which have associated images
    cal_focal_ratios : list(float), optional
        Output focal ratios which were used as calibration images
    """
    def __init__(self,
                 name,
                 folder,
                 new=True,
                 input_focal_ratios=[2.5, 3.0, 3.5, 4.0, 4.5, 5.0],
                 cal_focal_ratios=[3.0, 4.0, 5.0]):
        self.name = name
        self.new = new
        self.folder = folder
        self.input_focal_ratios = input_focal_ratios
        self.cal_focal_ratios = cal_focal_ratios

def coreExtensionFRD(test_data):
    """Collects all relevant FRD info for core extension fibers

    Args:
        test_data [Data]: object containing test information

    Returns:
        output_dict:
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
            magnification : tuple(float, float, list(float))
                the magnification, standard deviation of the magnification, and
                a list of each measured magnification for the far field camera
            name : string
                name given to the test
    """
    name = test_data.name
    new = test_data.new
    folder = test_data.folder
    input_focal_ratios = test_data.input_focal_ratios
    cal_focal_ratios = test_data.cal_focal_ratios

    calibration = Calibration([folder + 'Dark/im_' + str(i).zfill(3) + '.fit' for i in xrange(10)],
                              None,
                              [folder + 'Ambient/im_' + str(i).zfill(3) + '.fit' for i in xrange(10)])

    magn_list = []
    for f in cal_focal_ratios:  
        images = [folder + 'Output ' + str(f) + '/im_' + str(i).zfill(3) + '.fit' for i in xrange(10)]
        im_obj = ImageAnalysis(images, calibration, magnification=1, threshold=2000)
        diameter = im_obj.getFiberDiameter(method='edge', units='microns', show_image=False)
        magn_list.append(diameter / ((4.0 / f) * 25400))
        saveArray(im_obj.getImage(), folder + 'Output ' + str(f) + ' Image.tif')
    magnification = np.array(magn_list).mean()
    stdev = np.array(magn_list).std()
    print 'Magnification:', magnification, '+/-', stdev

    fig1, ax1 = plt.subplots()
    encircled_energy_focal_ratios = []
    encircled_energy = []
    energy_loss = []
    output_focal_ratios = []
    for f in input_focal_ratios:
        print 'F/' + str(f) + ':'
        images = [folder + 'Input ' + str(f) + '/im_' + str(i).zfill(3) + '.fit' for i in xrange(10)]
        ff_obj = ImageAnalysis(images, calibration, magnification=magnification)
        frd_output = FRD(ff_obj, input_focal_ratio=f, focal_lim=(2.3, 6.0), res=0.1)

        encircled_energy_focal_ratios.append(frd_output[0])
        encircled_energy.append(frd_output[1])
        ax1.plot(frd_output[0], frd_output[1], label=str(f))

        energy_loss.append(frd_output[2])
        print 'Energy loss:', frd_output[2], '%'

        output_focal_ratios.append(frd_output[3])
        print 'Output F/#:', frd_output[3]

        saveArray(ff_obj.getImage(), folder + 'Input ' + str(f) + ' Image.tif')
        saveCrossSections(ff_obj.getImage(),
                          ff_obj.getFiberCentroid()[0],
                          ff_obj.getFiberCentroid()[1],
                          folder + 'Input ' + str(f) + ' Cross Sections.tif')

    ax1.set_xlabel('Output f/#')
    ax1.set_ylabel('Encircled Energy')
    ax1.legend(title='Input f/#', loc=3)
    ax1.grid()
    ax1.set_title('FRD: ' + name)
    fig1.savefig(folder + name + ' FRD.png')

    print

    focal_ratio_relative_error = stdev / np.sqrt(len(magn_list)) / magnification

    output_dict = {'input_focal_ratios': np.array(input_focal_ratios),
                   'encircled_energy_focal_ratios': np.array(encircled_energy_focal_ratios),
                   'encircled_energy': encircled_energy,
                   'energy_loss': energy_loss,
                   'output_focal_ratios': output_focal_ratios,
                   'magnification': (magnification, stdev, magn_list),
                   'error': focal_ratio_relative_error,
                   'name': name}

    with open(folder + 'Info.txt', 'w') as file:
        file.write(str(output_dict))

    return output_dict

if __name__ == '__main__':
    if platform == 'darwin':
        base_folder = '/home/ryanp/Fiber_Characterization/'
    else:
        base_folder = 'C:/Libraries/Box Sync/ExoLab/Fiber_Characterization/Image Analysis/data/frd/'

    ref = Data('Reference Fiber',
               base_folder+'Core Extension/2016-08-10 Reference Octagonal/',
               False)
    prototype_1 = Data('Prototype Fiber 1',
                       base_folder+'Core Extension/2016-08-04 Reference Octagonal/',
                       False,
                       [3.0, 3.5, 4.0, 4.5, 5.0],
                       [3.5])
    prototype_2 = Data('Prototype Fiber 2',
                        base_folder+'Core Extension/2016-08-09 Prototype Core Extension 2/',
                        False)
    prototype_A2 = Data('Prototype Fiber A2',
                        base_folder+'Core Extension/2017-01-11 Prototype A2/',
                        False)
    prototype_A3 = Data('Prototype Fiber A3',
                        base_folder+'Core Extension/2017-01-12 Prototype A3/',
                        True)

    tests = [ref, prototype_1, prototype_2, prototype_A2, prototype_A3]

    for test in tests:
        print test.name + ':'
        if test.new:
            test.info = coreExtensionFRD(test)
        else:
            with open(test.folder + 'info.txt') as file:
                test.info = literal_eval(file.read())

    plt.figure()
    for test in tests:
        output = test.info
        plt.errorbar(output['input_focal_ratios'],
                     output['energy_loss'],
                     xerr=output['error']*output['input_focal_ratios'],
                     label=output['name'])
    plt.xlabel('Input f/#')
    plt.ylabel('Energy Loss [%]')
    plt.grid()
    plt.legend(loc=4)
    plt.title('Energy Loss at Constant F/#')
    plt.savefig(base_folder + 'Core Extension/Energy Loss.png')

    plt.figure()
    for test in tests:
        output = test.info
        plt.errorbar(output['input_focal_ratios'],
                     output['output_focal_ratios'],
                     xerr=output['error']*output['input_focal_ratios'],
                     yerr=output['error']*output['output_focal_ratios'],
                     label=output['name'])
    plt.plot(prototype2_2['input_focal_ratios'], prototype2_2['input_focal_ratios'], label='Ideal', linestyle='--', color='black')
    plt.xlabel('Input f/#')
    plt.ylabel('Output f/#')
    plt.grid()
    plt.legend(loc=2)
    plt.title('FRD Comparison')
    plt.savefig(base_folder + 'Core Extension/Input vs Output.png')

    plt.figure()
    for tests in tests:
        output = test.info
        for i, focal_ratio in enumerate([2.5, 3.0, 3.5, 4.0, 4.5, 5.0]):
            if focal_ratio in output['input_focal_ratios']:
                plt.subplot(3, 2, i+1)
                index = output['input_focal_ratios'].index(focal_ratio)
                plt.plot(output['encircled_energy_focal_ratios'][index], output['encircled_energy'][index], label=output['name'], linewidth=1)
                plt.xlabel('Output f/#', fontsize=8)
                plt.ylabel('Encircled Energy', fontsize=8)
                plt.xticks(fontsize=8)
                plt.yticks(fontsize=8)
                plt.grid()
                plt.legend(loc=3, fontsize=8)
                plt.title('Input f/# = ' + str(focal_ratio), fontsize=8)
    plt.savefig(base_folder + 'Core Extension/Encircled Energy Comparison.png')

