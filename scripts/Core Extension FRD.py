from ImageAnalysis import ImageAnalysis
from Calibration import Calibration
from FiberProperties import FRD
from NumpyArrayHandler import saveArray, saveCrossSections
import matplotlib.pyplot as plt
import numpy as np
from ast import literal_eval
from sys import platform

def coreExtensionFRD(folder, name, input_focal_ratios, test_focal_ratios):
    """Collects all relevant FRD info for core extension fibers

    Args:
        folder: top level folder where images are contained
        name: name to be used for saving plots

    Returns:
        output_dict:
            input_focal_ratio [list(float)]: list of the given input focal
                ratios
            encircled_energy [list(list(float))]: list of the lists of
                encircled energies for each input focal ratio
            encircled_energy_focal_ratios [list(list(float))]: independent
                variable (output f/#) corresponding to each encircled energy
            energy_loss [list(float)]: list of energy losses for each
                input focal ratio
            output_focal_ratios [list(float)]: list of calculated output focal
                ratio for each input focal ratio
            magnification [tuple(float, float, list(float))]: the magnification,
                standard deviation of the magnification, and a list of each
                measured magnification for the far field camera
            name [string]: name given to the test
    """
    calibration = Calibration([folder + 'Dark/im_' + str(i).zfill(3) + '.fit' for i in xrange(10)],
                              None,
                              [folder + 'Ambient/im_' + str(i).zfill(3) + '.fit' for i in xrange(10)])

    magn_list = []
    for f in test_focal_ratios:  
        images = [folder + 'Output ' + str(f) + '/im_' + str(i).zfill(3) + '.fit' for i in xrange(10)]
        im_obj = ImageAnalysis(images, calibration, magnification=1, threshold=10000)
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

    output_dict = {'input_focal_ratios': input_focal_ratios,
                   'encircled_energy_focal_ratios': encircled_energy_focal_ratios,
                   'encircled_energy': encircled_energy,
                   'energy_loss': energy_loss,
                   'output_focal_ratios': output_focal_ratios,
                   'magnification': (magnification, stdev, magn_list),
                   'name': name}

    with open(folder + 'Info.txt', 'w') as file:
        file.write(str(output_dict))

    return output_dict

if __name__ == '__main__':
    if platform == 'darwin':
        base_folder = '/home/ryanp/Fiber_Characterization/'
    else:
        base_folder = 'C:/Libraries/Box Sync/ExoLab/Fiber_Characterization/'
    ref_1 = dict(name='Reference Fiber Trial 1',
                 new_data=True,
                 folder=base_folder+'FRD Measurements/Core Extension/2016-08-04 Reference Octagonal/',
                 input_focal_ratios=[3.0, 3.5, 4.0, 4.5, 5.0],
                 cal_focal_ratios=[3.5])
    ref_2 = dict(name='Reference Fiber Trial 2',
                 new_data=True,
                 folder=base_folder+'FRD Measurements/Core Extension/2016-08-10 Reference Octagonal/',
                 input_focal_ratios=[2.5, 3.0, 3.5, 4.0, 4.5, 5.0],
                 cal_focal_ratios=[3.0, 4.0, 5.0])
    prototype1 = dict(name='Prototype Fiber 1 Trial 1',
                      new_data=True,
                      folder=base_folder+'FRD Measurements/Core Extension/2016-08-05 Prototype Core Extension 1/',
                      input_focal_ratios=[3.0, 3.5, 4.0, 4.5, 5.0],
                      cal_focal_ratios=[3.5])
    prototype2_1 = dict(name='Prototype Fiber 2 Trial 1',
                        new_data=True,
                        folder=base_folder+'FRD Measurements/Core Extension/2016-08-08 Prototype Core Extension 2/',
                        input_focal_ratios=[3.0, 3.5, 4.0, 4.5, 5.0],
                        cal_focal_ratios=[3.5, 5.0])
    prototype2_2 = dict(name='Prototype Fiber 2 Trial 2',
                        new_data=True,
                        folder=base_folder+'FRD Measurements/Core Extension/2016-08-09 Prototype Core Extension 2/',
                        input_focal_ratios=[2.5, 3.0, 3.5, 4.0, 4.5, 5.0],
                        cal_focal_ratios=[3.0, 4.0, 5.0])

    for test in [ref_1, ref_2, prototype1, prototype2_1, prototype2_2]:
        print test['name'] + ':'
        if test['new_data']:
            test['info'] = coreExtensionFRD(test['folder'],
                                            test['name'],
                                            test['input_focal_ratios'],
                                            test['cal_focal_ratios'])
        else:
            with open(test['folder'] + 'info.txt') as file:
                test['info'] = literal_eval(file.read())

    plt.figure()
    for output in [ref_1, ref_2, prototype1, prototype2_1, prototype2_2]:
        output = output['info']
        plt.plot(output['input_focal_ratios'], output['energy_loss'], label=output['name'])
    plt.xlabel('Input f/#')
    plt.ylabel('Energy Loss [%]')
    plt.grid()
    plt.legend(loc=4)
    plt.title('Energy Loss at Constant F/#')
    plt.savefig('FRD Measurements/Core Extension/Energy Loss.png')

    plt.figure()
    for output in [ref_1, ref_2, prototype1, prototype2_1, prototype2_2]:
        output = output['info']
        plt.plot(output['input_focal_ratios'], output['output_focal_ratios'], label=output['name'])
    plt.plot(prototype2_2['input_focal_ratios'], prototype2_2['input_focal_ratios'], label='Ideal', linestyle='--', color='black')
    plt.xlabel('Input f/#')
    plt.ylabel('Output f/#')
    plt.grid()
    plt.legend(loc=2)
    plt.title('FRD Comparison')
    plt.savefig('FRD Measurements/Core Extension/Input vs Output.png')

    plt.figure()
    for output in [ref_1, ref_2, prototype1, prototype2_1, prototype2_2]:
        output = output['info']
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
    plt.savefig('FRD Measurements/Core Extension/Encircled Energy Comparison.png')

